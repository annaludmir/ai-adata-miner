#!/usr/bin/env python3
"""22 - Gene programmes at single-cell level: per-cell scores and within-cell co-expression.

Pseudobulks say how much of a list a cell type expresses on average. They
cannot say whether every cell expresses it a little or a subset expresses it
a lot, nor whether the genes co-vary from cell to cell inside one cell type.
This pass measures both, in one read of X.

A. Per-cell programme scores (step 3 part B3 and B7)
   Score = mean log1p(CP10K) of the programme's genes minus the mean of
   expression-matched control genes (for each member, CTRL_PER_GENE genes from
   its bin of total UMIs; Tirosh et al. 2016 / Seurat AddModuleScore). Scored:
   user gene lists, seed NDD and cell-cycle panels, cell-state programmes
   (deep / upper layer, excitatory / inhibitory, OPC, astrocyte precursor,
   pre-OPC) and state contrasts (deep - upper, excitatory - inhibitory,
   OPC - astrocyte precursor). Each list and panel also gets N_RANDOM_PROGRAMMES
   random programmes drawing each member from the same expression bin, scored
   the same way -- the reference for "is the list's spread across cells unusual".
   Per group (cell class; class x age; class x phase) the pass keeps each score's
   count, mean, SD and a histogram on fixed edges.

B. Within-cell co-expression (part B4)
   Per cell class, up to CO_MAX_CELLS randomly chosen cells: covariance of
   log1p(CP10K) over a gene pool (list and panel genes, plus random expressed
   genes as the null pool) together with log total UMIs, then correlations
   with depth partialled out. For each list / panel: coherence = mean pairwise
   correlation of its genes (genes with mean log1p(CP10K) >= CO_MIN_MEAN in the
   class), each gene's connectivity, and a null from CO_N_NULL random sets of
   pool genes matched on their mean in that class. The random programmes of
   part A are tested the same way, so step 3 can check the null's calibration.

Needs script 09's gene_selection.csv (gene totals) for the same stratum.

Outputs (csv_exports/<dataset>__<chem>/22_cell_programs/)
  programs.csv                 every scored programme: kind, source, genes, controls
  program_score_summary.csv    grouping x group x programme: n_cells, mean, sd
  program_score_hist.csv       ... histogram counts on fixed edges (nonzero bins)
  within_cell_coherence.csv    class x list: coherence, null mean / sd, effect, p
  within_cell_connectivity.csv class x list x gene: mean correlation with the other members
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import scipy.sparse as sp

import config
from lib import cli
from lib.aggregate import combine_keys, drop_empty_levels, group_codes
from lib.io_utils import (Manifest, XReader, add_derived_obs_columns, chemistry_mask,
                          exclusion_mask, gene_frame, log, read_obs, read_var, resolve_role)
from lib.panels import load_panels

SCRIPT = "22_cell_level_programs"
SUBDIR = "22_cell_programs"
STATES = {
    "deep_layer": ["BCL11B", "TBR1", "FEZF2", "SOX5"],
    "upper_layer": ["SATB2", "CUX2", "POU3F2", "POU3F3"],
    "excitatory": ["SLC17A6", "SLC17A7", "NEUROD2", "NEUROD6"],
    "inhibitory": ["GAD1", "GAD2", "SLC32A1", "DLX5"],
    "opc": ["PDGFRA", "CSPG4", "OLIG1", "SOX10"],
    "astrocyte_precursor": ["AQP4", "GFAP", "S100B", "ALDH1L1", "SPARCL1"],
    "pre_opc": ["EGFR", "ASCL1", "OLIG2", "DLL3"],
}
CONTRASTS = {"deep_vs_upper": ("deep_layer", "upper_layer"),
             "excitatory_vs_inhibitory": ("excitatory", "inhibitory"),
             "opc_vs_astrocyte": ("opc", "astrocyte_precursor")}
PANEL_GROUPS = ("user_lists", "ndd", "cell_cycle")
MIN_PROGRAM_GENES = 3
N_EXPR_BINS = 25
CTRL_PER_GENE = 5
N_RANDOM_PROGRAMMES = 3
EDGES = np.round(np.arange(-2.0, 4.0001, 0.05), 4)    # histogram edges for scores (log1p units)
CO_MAX_CELLS = 20_000
CO_MAX_MEMBERS = 2_500
CO_NULL_POOL = 1_000
CO_MIN_MEAN = 0.02
CO_MIN_CELLS = 200
CO_N_NULL = 500
CO_MIN_GENES = 5


def map_members(genes: list[str], sym_upper: dict, ensg: dict) -> list[int]:
    out = []
    for g in genes:
        g = str(g).strip()
        i = sym_upper.get(g.upper(), ensg.get(g.split(".")[0]))
        if i is not None:
            out.append(i)
    return sorted(set(out))


def build_programmes(genes: pd.DataFrame, totals: np.ndarray, rng: np.random.Generator):
    """[(name, kind, source, member idx, control idx)] and the expression bins."""
    sym = genes["symbol"].astype(str).to_numpy()
    sym_upper = {}
    for i, s in enumerate(sym):
        sym_upper.setdefault(s.upper(), i)
    ensg = {}
    if "accession_base" in genes:
        for i, a in enumerate(genes["accession_base"].astype(str)):
            ensg.setdefault(a, i)
    expressed = totals > 0
    bins = np.full(len(totals), -1)
    lv = np.log1p(totals[expressed])
    bins[expressed] = np.digitize(lv, np.quantile(lv, np.linspace(0, 1, N_EXPR_BINS + 1)[1:-1]))
    pools = {b: np.nonzero(bins == b)[0] for b in range(N_EXPR_BINS)}

    def controls(members: list[int]) -> list[int]:
        mem = set(members)
        out = set()
        for m in members:
            pool = [x for x in pools[bins[m]] if x not in mem]
            if pool:
                out.update(rng.choice(pool, size=min(CTRL_PER_GENE, len(pool)), replace=False).tolist())
        return sorted(out)

    progs = []
    panels = load_panels()
    for grp in PANEL_GROUPS:
        for pname, glist in panels.get(grp, {}).items():
            mem = [i for i in map_members(glist, sym_upper, ensg) if expressed[i]]
            if len(mem) < MIN_PROGRAM_GENES:
                continue
            name = f"{'list' if grp == 'user_lists' else 'seed'}:{pname}"
            ctrl = controls(mem)
            progs.append((name, "programme", name, mem, ctrl))
            for r in range(N_RANDOM_PROGRAMMES):
                rmem = sorted({int(rng.choice(pools[bins[m]])) for m in mem})
                progs.append((f"{name}~random{r + 1}", "random", name, rmem, ctrl))
    for sname, glist in STATES.items():
        mem = [i for i in map_members(glist, sym_upper, ensg) if expressed[i]]
        if len(mem) >= 2:
            progs.append((f"state:{sname}", "state", sname, mem, controls(mem)))
    return progs, bins


def weight_matrix(progs, cols_pos: dict, contrasts: dict) -> tuple[sp.csr_matrix, list[str]]:
    rows, cols, vals, names = [], [], [], []
    by_name = {}
    for j, (name, kind, src, mem, ctrl) in enumerate(progs):
        by_name[name] = j
        for g in mem:
            rows.append(cols_pos[g]); cols.append(j); vals.append(1.0 / len(mem))
        for g in ctrl:
            rows.append(cols_pos[g]); cols.append(j); vals.append(-1.0 / len(ctrl))
        names.append(name)
    W = sp.coo_matrix((vals, (rows, cols)), shape=(len(cols_pos), len(progs))).tocsc()
    extra = []
    for cname, (a, b) in contrasts.items():
        if f"state:{a}" in by_name and f"state:{b}" in by_name:
            extra.append(W[:, by_name[f"state:{a}"]] - W[:, by_name[f"state:{b}"]])
            names.append(f"contrast:{cname}")
    if extra:
        W = sp.hstack([W] + extra).tocsc()
    return W.tocsr(), names


class ScoreAccumulator:
    def __init__(self, n_groups: int, n_scores: int):
        self.nb = len(EDGES) + 1
        self.n = np.zeros((n_groups, n_scores))
        self.s1 = np.zeros((n_groups, n_scores))
        self.s2 = np.zeros((n_groups, n_scores))
        self.hist = np.zeros(n_groups * n_scores * self.nb, dtype=np.int64)
        self.shape = (n_groups, n_scores)

    def update(self, S: np.ndarray, codes: np.ndarray) -> None:
        ok = codes >= 0
        if not ok.any():
            return
        S, c = S[ok], codes[ok]
        ng, ns = self.shape
        np.add.at(self.n, c, 1.0)
        np.add.at(self.s1, c, S)
        np.add.at(self.s2, c, S ** 2)
        b = np.searchsorted(EDGES, S, side="right")          # 0 .. len(EDGES)
        flat = (c[:, None] * ns + np.arange(ns)[None, :]) * self.nb + b
        self.hist += np.bincount(flat.ravel(), minlength=self.hist.size)


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)
    sel_path = config.CSV_EXPORTS / ns / "09_pseudobulk" / "gene_selection.csv"
    if not sel_path.exists():
        log("  09_pseudobulk/gene_selection.csv missing -- run 09 first; skipping")
        man.flush()
        return
    obs = add_derived_obs_columns(read_obs(path), key)
    keep = chemistry_mask(obs, chem)
    if keep is None:
        log(f"  no chemistry column -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    keep = keep.to_numpy() & exclusion_mask(obs, key, within=keep.to_numpy())[0]
    cls_col = resolve_role(obs, "cell_class")
    if cls_col is None:
        log("  no cell-class column; skipping")
        man.flush()
        return
    with XReader(path) as xr:
        n_obs = xr.n_obs
    n_cells = min(args.limit_cells, n_obs) if args.limit_cells else n_obs
    keep[n_cells:] = False
    genes = gene_frame(read_var(path), key)
    sel = pd.read_csv(sel_path)
    if len(sel) != len(genes):
        log("  gene_selection does not match this file's genes -- skipping")
        man.flush()
        return
    totals = sel["total_umis"].to_numpy(float)
    rng = np.random.default_rng(config.RANDOM_SEED)
    progs, bins = build_programmes(genes, totals, rng)
    log(f"  {sum(p[1] == 'programme' for p in progs)} lists / panels, "
        f"{sum(p[1] == 'state' for p in progs)} state programmes, "
        f"{sum(p[1] == 'random' for p in progs)} random programmes")

    # gene pool for within-cell co-expression: list / panel members + random expressed genes
    members = sorted({g for p in progs if p[1] == "programme" for g in p[3]},
                     key=lambda g: -totals[g])[:CO_MAX_MEMBERS]
    rand_members = sorted({g for p in progs if p[1] == "random" for g in p[3]} - set(members),
                          key=lambda g: -totals[g])[:CO_MAX_MEMBERS // 2]
    members = members + rand_members
    expressed = np.nonzero(totals > 0)[0]
    others = np.setdiff1d(expressed, members)
    null_pool = rng.choice(others, size=min(CO_NULL_POOL, others.size), replace=False).tolist()
    pool = sorted(set(members) | set(null_pool))
    needed = sorted({g for p in progs for g in p[3] + p[4]} | set(pool))
    cols_pos = {g: i for i, g in enumerate(needed)}
    W, score_names = weight_matrix(progs, cols_pos, CONTRASTS)
    pool_pos = np.array([cols_pos[g] for g in pool])
    log(f"  {len(needed):,} genes read per cell; co-expression pool {len(pool):,} genes "
        f"({len(members) - len(rand_members):,} list / panel, {len(rand_members):,} random-programme, "
        f"{len(null_pool):,} null)")

    cls_series = obs[cls_col].astype("object")
    groupings = {"cell_class": cls_series}
    if "age_pcw" in obs:
        groupings["cell_class_x_age"] = combine_keys(obs.assign(_a=obs["age_pcw"].astype(str)), [cls_col, "_a"])
    if "cyclephase_h" in obs and obs["cyclephase_h"].notna().any():
        groupings["cell_class_x_phase"] = combine_keys(obs, [cls_col, "cyclephase_h"])
    accs, codes_map, levels_map = {}, {}, {}
    for name, series in groupings.items():
        codes, levels = group_codes(series)
        codes = np.where(keep, codes, -1)
        codes, levels = drop_empty_levels(codes, levels)
        codes_map[name], levels_map[name] = codes, levels
        accs[name] = ScoreAccumulator(len(levels), len(score_names))

    # co-expression cell sample per class
    ccodes, clevels = codes_map["cell_class"], levels_map["cell_class"]
    co_pick = np.zeros(len(obs), dtype=bool)
    for k in range(len(clevels)):
        idx = np.nonzero(ccodes == k)[0]
        if idx.size >= CO_MIN_CELLS:
            co_pick[rng.choice(idx, size=min(CO_MAX_CELLS, idx.size), replace=False)] = True
    G = len(pool)
    co = {k: {"n": 0.0, "sx": np.zeros(G), "sxx": np.zeros((G, G)), "su": 0.0, "suu": 0.0,
              "sxu": np.zeros(G)} for k in range(len(clevels)) if (co_pick & (ccodes == k)).any()}
    log(f"  co-expression on {int(co_pick.sum()):,} cells in {len(co)} classes")

    with XReader(path) as xr:
        for start, stop, chunk in xr.iter_chunks(args.chunk_size):
            if start >= n_cells:
                break
            stop = min(stop, n_cells)
            chunk = chunk[: stop - start]
            m = keep[start:stop]
            if not m.any():
                continue
            tot = np.asarray(chunk.sum(axis=1)).ravel()
            scale = np.divide(config.TARGET_SUM, tot, out=np.zeros_like(tot, dtype=float), where=tot > 0)
            L = (sp.diags(scale) @ chunk[:, needed]).tocsr()
            L.data = np.log1p(L.data)
            S = np.asarray((L @ W).todense()) if sp.issparse(L @ W) else L @ W
            for name, acc in accs.items():
                acc.update(S, codes_map[name][start:stop])
            pick = co_pick[start:stop]
            if pick.any():
                Lp = L[pick][:, pool_pos]
                u = np.log1p(tot[pick])
                cc = ccodes[start:stop][pick]
                for k in np.unique(cc):
                    if k not in co:
                        continue
                    r = cc == k
                    M = Lp[r]
                    d = co[k]
                    d["n"] += r.sum()
                    d["sx"] += np.asarray(M.sum(axis=0)).ravel()
                    d["sxx"] += np.asarray((M.T @ M).todense())
                    d["su"] += u[r].sum()
                    d["suu"] += (u[r] ** 2).sum()
                    d["sxu"] += np.asarray(M.T @ u[r]).ravel()
            log(f"    {stop:,} / {n_cells:,} cells")

    # ---- programmes and score summaries ------------------------------------------
    sym = genes["symbol"].astype(str).to_numpy()
    pt = pd.DataFrame([{"program": p[0], "kind": p[1], "source": p[2], "n_genes": len(p[3]),
                        "n_controls": len(p[4]), "genes": "|".join(sym[p[3]][:200])} for p in progs])
    pt = pd.concat([pt, pd.DataFrame([{"program": n, "kind": "contrast", "source": n.split(":")[1],
                                       "n_genes": np.nan, "n_controls": np.nan, "genes": ""}
                                      for n in score_names if n.startswith("contrast:")])], ignore_index=True)
    man.write(pt, "programs", "Scored programmes: kind (programme / random / state / contrast), genes, controls",
              subdir=SUBDIR)
    summ, hist = [], []
    lo = np.concatenate([[-np.inf], EDGES])
    for name, acc in accs.items():
        levels = levels_map[name]
        ng, nsc = acc.shape
        H = acc.hist.reshape(ng, nsc, acc.nb)
        for g in range(ng):
            n = acc.n[g]
            if n.max() == 0:
                continue
            mean = acc.s1[g] / np.maximum(n, 1)
            sd = np.sqrt(np.maximum(acc.s2[g] / np.maximum(n, 1) - mean ** 2, 0))
            summ.append(pd.DataFrame({"grouping": name, "group": levels[g], "program": score_names,
                                      "n_cells": n.astype(int), "mean": mean, "sd": sd}))
            j, b = np.nonzero(H[g])
            hist.append(pd.DataFrame({"grouping": name, "group": levels[g],
                                      "program": np.array(score_names)[j], "bin_lo": lo[b],
                                      "count": H[g][j, b]}))
    man.write(pd.concat(summ, ignore_index=True), "program_score_summary",
              "Per grouping x group x programme: cells, mean and SD of the per-cell score", subdir=SUBDIR)
    man.write(pd.concat(hist, ignore_index=True), "program_score_hist",
              f"Per grouping x group x programme: cells per score bin (left edge; width {EDGES[1] - EDGES[0]:g}; "
              "-inf = below the first edge)", subdir=SUBDIR)

    # ---- within-cell co-expression ------------------------------------------------
    coh_rows, con_rows = [], []
    pool_arr = np.array(pool)
    pidx = {g: i for i, g in enumerate(pool)}
    for k, d in co.items():
        n = d["n"]
        if n < CO_MIN_CELLS:
            continue
        mx = d["sx"] / n
        cov = d["sxx"] / n - np.outer(mx, mx)
        mu = d["su"] / n
        vu = d["suu"] / n - mu ** 2
        cxu = d["sxu"] / n - mx * mu
        cov_p = cov - np.outer(cxu, cxu) / vu if vu > 0 else cov     # depth partialled out
        sdp = np.sqrt(np.clip(np.diag(cov_p), 0, None))
        ok_g = (mx >= CO_MIN_MEAN) & (sdp > 0)
        R = np.divide(cov_p, np.outer(sdp, sdp), out=np.zeros_like(cov_p), where=np.outer(sdp, sdp) > 0)
        np.fill_diagonal(R, 1.0)
        okidx = np.nonzero(ok_g)[0]
        mbins = np.full(G, -1)
        mbins[okidx] = np.digitize(np.log1p(mx[okidx]),
                                   np.quantile(np.log1p(mx[okidx]), np.linspace(0, 1, 11)[1:-1]))
        mpools = {b: okidx[mbins[okidx] == b] for b in range(10)}
        for name, kind, src, mem, ctrl in progs:
            if kind not in ("programme", "random"):      # random programmes calibrate the null
                continue
            idx = np.array([pidx[g] for g in mem if g in pidx and ok_g[pidx[g]]])
            if idx.size < CO_MIN_GENES:
                continue
            kk = idx.size
            sub = R[np.ix_(idx, idx)]
            coh = float((sub.sum() - kk) / (kk * (kk - 1)))
            conn = (sub.sum(axis=1) - 1) / (kk - 1)
            B = np.zeros((CO_N_NULL, G))
            for s in range(CO_N_NULL):
                for b in np.unique(mbins[idx]):
                    need = int((mbins[idx] == b).sum())
                    pool_b = mpools[b]
                    pick = rng.choice(pool_b, size=min(need, pool_b.size), replace=False)
                    B[s, pick] = 1
            ks = B.sum(axis=1)
            null = ((B @ R) * B).sum(axis=1)
            null = (null - ks) / (ks * (ks - 1))
            mu0, sd0 = float(null.mean()), float(null.std())
            p = (np.sum(np.abs(null - mu0) >= abs(coh - mu0)) + 1) / (CO_N_NULL + 1)
            coh_rows.append({"cell_class": clevels[k], "program": name, "kind": kind, "n_genes": kk,
                             "n_cells": int(n),
                             "coherence": coh, "null_mean": mu0, "null_sd": sd0,
                             "effect_vs_null_sd": (coh - mu0) / sd0 if sd0 > 0 else np.nan, "perm_p": p})
            if kind != "programme":
                continue
            con_rows.append(pd.DataFrame({"cell_class": clevels[k], "program": name,
                                          "gene": sym[pool_arr[idx]], "connectivity": conn,
                                          "mean_log1p_cp10k": mx[idx]}))
    if coh_rows:
        man.write(pd.DataFrame(coh_rows), "within_cell_coherence",
                  "Per class x list: mean pairwise correlation across single cells (depth partialled), "
                  f"vs {CO_N_NULL} random pool sets matched on class mean", subdir=SUBDIR)
        man.write(pd.concat(con_rows, ignore_index=True), "within_cell_connectivity",
                  "Per class x list x gene: mean within-cell correlation with the other members", subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
