#!/usr/bin/env python3
"""24 - Cell subsample for Spectra (knowledge-guided gene programmes), one per stratum.

Spectra (Kunes et al. 2023) fits gene programmes to single cells, guided by
prior gene sets. It needs a cells x genes matrix in memory, so each stratum
gets a stratified subsample: up to N_PER_CLASS cells per cell class (classes
with >= MIN_CLASS_CELLS cells), at most MAX_CELLS in all, drawn with a fixed
seed so a rerun reproduces it.

Per subsample: log1p(CP10K) in .X, raw counts in .layers['counts'];
.obs = cell class, age, donor, region, cell-cycle phase, UMIs; .var indexed by
gene symbol (made unique), with 'highly_variable' = the N_HVG most dispersed
genes (scanpy, seurat flavour) -- Spectra uses those plus every gene in a
prior set (script 25).

Writes the subsample to <AIM_WORK>/spectra/<dataset>__<chem>/input.h5ad (not a
CSV export: it is a working file), and a summary CSV to
csv_exports/<dataset>__<chem>/24_spectra/input_summary.csv.
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
from lib.io_utils import (Manifest, XReader, add_derived_obs_columns, chemistry_mask, exclusion_mask,
                          gene_frame, log, read_obs, read_var, resolve_role)

SCRIPT = "24_spectra_input"
SUBDIR = "24_spectra"
N_PER_CLASS = 2500
MAX_CELLS = 25_000
MIN_CLASS_CELLS = 200
N_HVG = 3000


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    import anndata as ad
    import scanpy as sc
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)
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
    cls_col = resolve_role(obs, "cell_class")
    cls = obs[cls_col].astype(str).to_numpy()
    rng = np.random.default_rng(config.RANDOM_SEED)
    counts = pd.Series(cls[keep]).value_counts()
    classes = [c for c, k in counts.items() if k >= MIN_CLASS_CELLS]
    want = {c: min(N_PER_CLASS, int(counts[c])) for c in classes}
    total = sum(want.values())
    if total > MAX_CELLS:
        f = MAX_CELLS / total
        want = {c: max(MIN_CLASS_CELLS, int(k * f)) for c, k in want.items()}
    pick = []
    for c in classes:
        idx = np.nonzero(keep & (cls == c))[0]
        pick.append(rng.choice(idx, size=want[c], replace=False))
    pick = np.sort(np.concatenate(pick)) if pick else np.array([], dtype=int)
    log(f"  subsample: {len(pick):,} cells from {len(classes)} classes "
        + ", ".join(f"{c} {want[c]}" for c in classes))
    if len(pick) == 0:
        man.flush()
        return

    rows = []
    with XReader(path) as xr:
        for start, stop, chunk in xr.iter_chunks(args.chunk_size):
            if start >= n_cells:
                break
            sel = pick[(pick >= start) & (pick < stop)] - start
            if sel.size:
                rows.append(chunk[sel])
    X = sp.vstack(rows).tocsr()
    genes = gene_frame(read_var(path), key)
    sym = pd.Index(genes["symbol"].astype(str).to_numpy())
    sym = ad.utils.make_index_unique(sym)
    o = obs.iloc[pick]
    meta = pd.DataFrame({"cell_class": cls[pick]}, index=o.index.astype(str))
    for role, name in (("donor", "donor"), ("region", "region")):
        col = resolve_role(obs, role)
        if col:
            meta[name] = o[col].astype(str).to_numpy()
    for col in ("age_pcw", "cyclephase_h"):
        if col in o:
            meta[col] = o[col].to_numpy()
    meta["total_umis"] = np.asarray(X.sum(axis=1)).ravel()
    var = pd.DataFrame({"accession": genes["accession"].astype(str).to_numpy()}, index=sym)
    a = ad.AnnData(X=X.astype(np.float32), obs=meta, var=var)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=config.TARGET_SUM)
    sc.pp.log1p(a)
    detected = np.asarray((a.layers["counts"] > 0).sum(axis=0)).ravel()
    a = a[:, detected >= 10].copy()
    sc.pp.highly_variable_genes(a, n_top_genes=min(N_HVG, a.n_vars - 1), flavor="seurat")
    out_dir = config.WORK_DIR / "spectra" / ns
    out_dir.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out_dir / "input.h5ad", compression="gzip")
    log(f"  wrote {out_dir / 'input.h5ad'}: {a.n_obs:,} cells x {a.n_vars:,} genes "
        f"({int(a.var.highly_variable.sum())} highly variable)")
    summ = meta.groupby("cell_class").size().rename("n_cells").reset_index()
    summ["n_genes_detected_10_cells"] = a.n_vars
    summ["n_highly_variable"] = int(a.var.highly_variable.sum())
    summ["input_h5ad"] = str(out_dir / "input.h5ad")
    man.write(summ, "input_summary", "Spectra subsample: cells per class, genes kept, highly variable genes",
              subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
