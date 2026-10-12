#!/usr/bin/env python3
"""28 - Score every cell with the gene programmes fitted on the subsample.

The programme fits (Spectra, script 25; cNMF, 26; scHPF, 27) ran on a
stratified subsample of ~12,500-25,000 cells per stratum. Here each fit's gene
side is held fixed and every cell of the stratum -- all cells that pass the
same chemistry and QC filters -- gets its programme activities, so programme
activity over age, region and donor rests on all cells, not the subsample.

Per fit, the cell side is re-estimated the way the method itself does it:
  Spectra   the exact model from model.pt (programmes x gene scaling, global
            plus the cell's own class); loadings by Poisson maximum likelihood
            on log1p(CP10K), Spectra's data term, solved by multiplicative
            updates (N_ITER_SPECTRA); cell score = loading x the programme's
            mean weight, as Spectra reports it. Cells of a class the fit did not
            model get the global programmes only.
  cNMF      cNMF's own final step: usages refitted (its refit_usage, its NMF
            settings) on variance-scaled TPM of its genes against its TPM
            spectra, with the subsample's per-gene scale; usages sum to 1.
  scHPF     the package's project(): cell distributions fitted with the gene
            distributions fixed, in chunks.

Check: the subsample cells are among the projected ones, so their projected
activity is compared with the fitted activity programme by programme
(projection_check.csv, Pearson r; it should be close to 1).

Reads the stratum's h5ad (streamed), the fits' work files under <AIM_WORK>
and their CSV exports. Writes csv_exports/<ns>/28_projection/<fit>/
  factor_scores_by_group.csv   mean activity per programme x (class | class x age |
                               class x region x age | donor), over all cells
  projection_check.csv         per programme: r between projected and fitted activity
  projection_summary.csv       cells projected, cells per class, fit settings
<fit> = spectra_full, spectra_core, spectra_string, cnmf, schpf (those that exist).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import (Manifest, XReader, add_derived_obs_columns, chemistry_mask, exclusion_mask,
                          gene_frame, log, read_obs, read_var, resolve_role)

SCRIPT = "28_project_programs"
SUBDIR = "28_projection"
SPECTRA_FITS = {"spectra_full": ("25_spectra", ""), "spectra_core": ("25_spectra_core", "_core"),
                "spectra_string": ("25_spectra_string", "_string")}
SPECTRA_DELTA = 0.001        # est_spectra's default delta (script 25 leaves it)
N_ITER_SPECTRA = 150
SCHPF_CHUNK = 20_000


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def _softmax_rows(x):
    e = np.exp(x - x.max(axis=1, keepdims=True))
    return e / e.sum(axis=1, keepdims=True)


class SpectraProjector:
    """Spectra's fitted programmes (exact, from model.pt) and a Poisson loading fit for new cells."""

    def __init__(self, ns: str, sub: str, suffix: str, gene_pos: dict):
        import torch
        exp = config.CSV_EXPORTS / ns / sub
        work = config.WORK_DIR / "spectra" / ns
        self.ok = (exp / "factor_info.csv").exists() and (work / f"model{suffix}.pt").exists()
        if not self.ok:
            return
        info = pd.read_csv(exp / "factor_info.csv")
        vocab = list(pd.read_csv(exp / "factor_gene_weights.csv", nrows=0).columns[1:])
        sd = torch.load(work / f"model{suffix}.pt", map_location="cpu", weights_only=False)
        self.factors = list(info.factor)
        self.cols = np.array([gene_pos.get(g, -1) for g in vocab])
        self.missing = int((self.cols < 0).sum())
        scopes = list(info.scope)
        # per scope: the scaled programmes (L x p), as Spectra builds `scaled` in __store_parameters
        self.S = {}
        for sc in dict.fromkeys(scopes):
            th = _softmax_rows(sd[f"theta.{sc}"].numpy().astype(np.float64))        # p x L, softmax over programmes
            gs = _sigmoid(sd[f"gene_scaling.{sc}"].numpy().astype(np.float64)) + SPECTRA_DELTA
            self.S[sc] = (th * gs[:, None]).T                                       # L x p
        self.idx = {sc: [i for i, s in enumerate(scopes) if s == sc] for sc in self.S}
        self.mean_w = np.zeros(len(scopes))
        for sc in self.S:
            self.mean_w[self.idx[sc]] = self.S[sc].mean(axis=1)
        self.fitted = None
        npz = work / f"cell_scores{suffix}.npz"
        if npz.exists():
            z = np.load(npz, allow_pickle=True)
            self.fitted = pd.DataFrame(z["scores"], index=pd.Index(z["cells"].astype(str)), columns=self.factors)

    def project(self, counts, totals: np.ndarray, classes: np.ndarray) -> np.ndarray:
        """counts: cells x all genes (csr); returns cells x factors cell scores."""
        X = np.zeros((counts.shape[0], len(self.cols)), dtype=np.float64)
        ok = self.cols >= 0
        X[:, ok] = counts[:, self.cols[ok]].toarray()
        X = np.log1p(X / np.maximum(totals, 1)[:, None] * config.TARGET_SUM)
        out = np.zeros((X.shape[0], len(self.factors)))
        for cls in np.unique(classes):
            rows = classes == cls
            scopes = ["global"] + ([cls] if cls in self.S else [])
            S = np.vstack([self.S[s] for s in scopes])
            fidx = sum((self.idx[s] for s in scopes), [])
            out[np.ix_(rows, fidx)] = self._loadings(X[rows], S) * self.mean_w[fidx]
        return out

    @staticmethod
    def _loadings(X: np.ndarray, S: np.ndarray) -> np.ndarray:
        """argmax_A sum(X log(AS) - AS), A >= 0: Poisson / KL multiplicative updates."""
        denom = S.sum(axis=1)                                            # K
        A = np.outer(X.sum(axis=1) / max(S.sum(), 1e-12), np.ones(S.shape[0])) + 1e-6
        for _ in range(N_ITER_SPECTRA):
            R = A @ S
            A *= ((X / np.maximum(R, 1e-12)) @ S.T) / denom
        return A


class CnmfProjector:
    """cNMF's last step (usage refit on variance-scaled TPM) for new cells."""

    def __init__(self, ns: str, gene_pos: dict):
        exp = config.CSV_EXPORTS / ns / "26_cnmf"
        work = config.WORK_DIR / "cnmf" / ns
        self.ok = (exp / "fit_summary.csv").exists() and (work / "cnmf").is_dir()
        if not self.ok:
            return
        import anndata as ad
        from cnmf import cNMF
        from cnmf.cnmf import load_df_from_npz
        k = int(pd.read_csv(exp / "fit_summary.csv").k.iloc[0])
        dt = str(pd.read_csv(exp / "fit_summary.csv").density_threshold.iloc[0]).replace(".", "_")
        self.obj = cNMF(output_dir=str(work), name="cnmf")
        tpm_genes = list(ad.read_h5ad(work / "counts.h5ad", backed="r").var_names)
        hvgs = [g for g in open(self.obj.paths["nmf_genes_list"]).read().split("\n") if g]
        std = load_df_from_npz(self.obj.paths["tpm_stats"])["__std"]
        spectra_tpm = load_df_from_npz(self.obj.paths["gene_spectra_tpm"] % (k, dt))
        self.H = spectra_tpm.loc[:, hvgs].div(std.loc[hvgs], axis=1)
        self.std = std.loc[hvgs].to_numpy()
        self.tpm_cols = np.array([gene_pos[g] for g in tpm_genes if g in gene_pos])
        self.hvg_cols = np.array([gene_pos[g] for g in hvgs])
        self.factors = [f"P{i:02d}" for i in self.H.index.astype(int)]
        self.H.index = self.factors
        usage, _, _, _ = self.obj.load_results(K=k, density_threshold=float(dt.replace("_", ".")))
        usage.columns = self.factors
        self.fitted = usage

    def project(self, counts, totals, classes) -> np.ndarray:
        tot = np.asarray(counts[:, self.tpm_cols].sum(axis=1)).ravel()
        X = counts[:, self.hvg_cols].toarray() / np.maximum(tot, 1)[:, None] * 1e6 / self.std
        U = np.asarray(self.obj.refit_usage(X, self.H.astype(X.dtype)))
        s = U.sum(axis=1, keepdims=True)
        return np.divide(U, s, out=np.zeros_like(U), where=s > 0)


class SchpfProjector:
    """scHPF's project(): cell distributions with the gene distributions fixed."""

    def __init__(self, ns: str, gene_pos: dict):
        exp = config.CSV_EXPORTS / ns / "27_schpf"
        work = config.WORK_DIR / "schpf" / ns
        self.ok = (exp / "factor_info.csv").exists() and (work / "model.joblib").exists()
        if not self.ok:
            return
        import joblib
        self.model = joblib.load(work / "model.joblib")
        genes = list(pd.read_csv(exp / "factor_gene_weights.csv", nrows=0).columns[1:])
        self.cols = np.array([gene_pos[g] for g in genes])
        self.factors = list(pd.read_csv(exp / "factor_info.csv").factor)
        src = config.WORK_DIR / "spectra" / ns / "input.h5ad"
        self.fitted = None
        if src.exists():
            import anndata as ad
            a = ad.read_h5ad(src, backed="r")
            X = a[:, genes].to_memory().layers["counts"]
            keep = np.asarray(X.sum(axis=1)).ravel() > 0           # the cells script 27 fitted
            self.fitted = pd.DataFrame(self.model.cell_score(), index=a.obs_names[keep].astype(str),
                                       columns=self.factors)

    def project(self, counts, totals, classes) -> np.ndarray:
        import scipy.sparse as sp
        X = sp.csr_matrix(counts[:, self.cols]).astype(np.float64)
        out = np.zeros((X.shape[0], len(self.factors)))
        has = np.asarray(X.sum(axis=1)).ravel() > 0
        rows = np.nonzero(has)[0]
        for i in range(0, rows.size, SCHPF_CHUNK):
            r = rows[i:i + SCHPF_CHUNK]
            proj = self.model.project(X[r].tocoo(), replace=False)
            out[r] = proj.cell_score()
        return out


class GroupSums:
    """Running per-group sums of activity, for class / class x age / class x region x age / donor."""

    def __init__(self, factors: list[str]):
        self.factors, self.sums, self.n = factors, {}, {}

    def add(self, meta: pd.DataFrame, scores: np.ndarray) -> None:
        cls = meta.cell_class.astype(str)
        groups = {"cell_class": cls}
        if "age_pcw" in meta:
            groups["cell_class_x_age"] = cls + " | " + meta.age_pcw.astype(str)
            if "region" in meta:
                groups["cell_class_x_region_x_age"] = (cls + " | " + meta.region.astype(str) + " | "
                                                       + meta.age_pcw.astype(str))
        if "donor" in meta:
            groups["donor"] = meta.donor.astype(str)
        S = pd.DataFrame(scores, columns=self.factors)
        for g, lab in groups.items():
            s = S.groupby(lab.to_numpy()).sum()
            n = lab.value_counts()
            self.sums[g] = s if g not in self.sums else self.sums[g].add(s, fill_value=0)
            self.n[g] = n if g not in self.n else self.n[g].add(n, fill_value=0)

    def table(self) -> pd.DataFrame:
        long = []
        for g, s in self.sums.items():
            m = s.div(self.n[g].reindex(s.index), axis=0)
            t = m.reset_index(names="group").melt(id_vars="group", var_name="factor", value_name="mean_score")
            t.insert(0, "grouping", g)
            t["n_cells"] = t.group.map(self.n[g]).astype(int)
            long.append(t)
        return pd.concat(long, ignore_index=True)


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    import anndata as ad
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    man = Manifest(ns, SCRIPT)
    path = cli.resolve_h5ad(key)
    obs = add_derived_obs_columns(read_obs(path), key)
    keep = chemistry_mask(obs, chem)
    if keep is None:
        log(f"  no chemistry column -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    keep = keep.to_numpy() & exclusion_mask(obs, key, within=keep.to_numpy())[0]
    with XReader(path) as xr:
        n_obs = xr.n_obs
    n_cells = min(args.limit_cells, n_obs) if args.limit_cells else n_obs
    keep[n_cells:] = False
    genes = gene_frame(read_var(path), key)
    sym = ad.utils.make_index_unique(pd.Index(genes["symbol"].astype(str).to_numpy()))
    gene_pos = {g: i for i, g in enumerate(sym)}
    meta = pd.DataFrame({"cell_class": obs[resolve_role(obs, "cell_class")].astype(str).to_numpy()},
                        index=obs.index.astype(str))
    for role in ("donor", "region"):
        col = resolve_role(obs, role)
        if col:
            meta[role] = obs[col].astype(str).to_numpy()
    if "age_pcw" in obs:
        meta["age_pcw"] = obs["age_pcw"].to_numpy()

    projectors = {name: SpectraProjector(ns, sub, suf, gene_pos) for name, (sub, suf) in SPECTRA_FITS.items()}
    projectors["cnmf"] = CnmfProjector(ns, gene_pos)
    projectors["schpf"] = SchpfProjector(ns, gene_pos)
    projectors = {k: p for k, p in projectors.items() if p.ok and (args.fits is None or k in args.fits)}
    if not projectors:
        log("  no finished programme fits for this stratum; skipping")
        man.flush()
        return
    log(f"  projecting {int(keep.sum()):,} cells onto: " + ", ".join(
        f"{k} ({len(p.factors)} programmes)" for k, p in projectors.items()))
    sums = {k: GroupSums(p.factors) for k, p in projectors.items()}
    check_rows = {k: [] for k in projectors}
    check_scores = {k: [] for k in projectors}
    targets = {k: set(p.fitted.index) if p.fitted is not None else set() for k, p in projectors.items()}
    done = 0
    with XReader(path) as xr:
        for start, stop, chunk in xr.iter_chunks(args.chunk_size):
            if start >= n_cells:
                break
            sel = np.nonzero(keep[start:stop])[0]
            if not sel.size:
                continue
            counts = chunk[sel].tocsr()
            totals = np.asarray(counts.sum(axis=1)).ravel()
            m = meta.iloc[start + sel]
            for k, p in projectors.items():
                sc = p.project(counts, totals, m.cell_class.to_numpy())
                sums[k].add(m.reset_index(drop=True), sc)
                hit = m.index.isin(list(targets[k]))
                if hit.any():
                    check_rows[k].extend(m.index[hit])
                    check_scores[k].append(sc[hit])
            done += sel.size
            log(f"  {done:,} / {int(keep.sum()):,} cells")
    for k, p in projectors.items():
        sub = f"{SUBDIR}/{k}"
        man.write(sums[k].table(), "factor_scores_by_group", "Mean projected activity per programme x class, class "
                  "x age, class x region x age and donor, over all cells", subdir=sub)
        chk = []
        if check_scores[k]:
            proj = pd.DataFrame(np.vstack(check_scores[k]), index=check_rows[k], columns=p.factors)
            fit = p.fitted.loc[proj.index]
            for f in p.factors:
                a, b = proj[f].to_numpy(), fit[f].to_numpy()
                ok = (a.std() > 0) and (b.std() > 0)
                chk.append({"factor": f, "n_cells": len(a), "pearson_r": float(np.corrcoef(a, b)[0, 1]) if ok else np.nan,
                            "mean_projected": float(a.mean()), "mean_fitted": float(b.mean())})
        chk = pd.DataFrame(chk)
        if len(chk):
            log(f"  {k}: projected vs fitted on {chk.n_cells.iloc[0]:,} subsample cells, median r "
                f"{chk.pearson_r.median():.3f} (min {chk.pearson_r.min():.3f})")
        man.write(chk, "projection_check", "Per programme: projected vs fitted activity on the subsample cells",
                  subdir=sub)
        per_class = sums[k].n["cell_class"]
        man.write(pd.DataFrame([{"fit": k, "cells_projected": int(per_class.sum()), "programmes": len(p.factors),
                                 "classes": "; ".join(f"{c} {int(v)}" for c, v in per_class.items()),
                                 "genes_missing_from_h5ad": int(getattr(p, "missing", 0))}]),
                  "projection_summary", "Cells projected and fit settings", subdir=sub)
    man.flush()


def main() -> None:
    p = cli.build_parser(__doc__)
    p.add_argument("--fits", nargs="*", default=None,
                   help="only these fits (spectra_full spectra_core spectra_string cnmf schpf)")
    args = p.parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
