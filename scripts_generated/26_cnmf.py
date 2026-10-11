#!/usr/bin/env python3
"""26 - Fit cNMF: gene programmes with no prior knowledge, per stratum.

Consensus NMF (Kotliar et al., eLife 2019; package cnmf) runs NMF many times
from random starts for each candidate number of programmes K, clusters the
replicate programmes and keeps their consensus. No gene set steers it, so it is
the prior-free counterpart of the Spectra fits (script 25): which programmes
appear without prior knowledge, and do the NDD lists gather in any of them?

Runs on the same cells and genes as Spectra: the stratified subsample of script
24 (<AIM_WORK>/spectra/<ns>/input.h5ad, raw counts) and its highly variable
genes. N_ITER replicates per K in K_GRID (CNMF_K_GRID="15 20 ..." overrides), spread over the job's CPUs. K is the
largest K whose stability (silhouette of the replicate programmes) is within
STABILITY_SLACK of the best -- the most resolved solution that is still about as
stable as any -- unless CNMF_K sets it; k_selection.csv keeps the statistics, so
another K needs only the consensus step (--k K --skip-factorize).

Outputs
  <AIM_WORK>/cnmf/<ns>/            cNMF working files (replicates, consensus)
  csv_exports/<ns>/26_cnmf/
    k_selection.csv             per K: stability, prediction error; chosen K
    factor_info.csv             per programme: top genes, nearest GO BP term, mean usage
                                (scope 'global', label 'new', free: the format of 25_spectra)
    factor_gene_weights.csv     programmes x genes: cNMF gene-spectra scores (z-score coefficients)
    factor_scores_by_group.csv  mean cell usage per programme x (class | class x age | class x region x age | donor)
    fit_summary.csv             K, replicates, cells, genes
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import Manifest, log
from lib.programs import write_group_scores

SCRIPT = "26_cnmf"
SUBDIR = "26_cnmf"
K_GRID = [15, 20, 25, 30, 40, 50, 60]
N_ITER = 100
DENSITY_THRESHOLD = 0.1       # drop replicate programmes far from their neighbours (cNMF's usual choice)
STABILITY_SLACK = 0.05
N_TOP = 50


THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def _worker(args):
    out_dir, name, i, n = args
    from cnmf import cNMF
    cNMF(output_dir=out_dir, name=name).factorize(worker_i=i, total_workers=n)


def choose_k(stats: pd.DataFrame) -> int:
    best = stats.silhouette.max()
    return int(stats.loc[stats.silhouette >= best - STABILITY_SLACK, "k"].max())


def go_terms() -> dict[str, set[str]]:
    gmt = Path(config.ANNOTATIONS_DIR) / "gmt" / "GO_Biological_Process_2023.gmt"
    if not gmt.exists():
        return {}
    out = {}
    for line in gmt.read_text().splitlines():
        parts = line.split("\t")
        genes = {g.upper() for g in parts[2:] if g}
        if 10 <= len(genes) <= 300:
            out[parts[0].split(" (GO:")[0]] = genes
    return out


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    import anndata as ad
    from multiprocessing import get_context
    from cnmf import cNMF
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    man = Manifest(ns, SCRIPT)
    src = config.WORK_DIR / "spectra" / ns / "input.h5ad"
    if not src.exists():
        log(f"  {src} missing -- run script 24 first; skipping")
        man.flush()
        return
    work = config.WORK_DIR / "cnmf" / ns
    work.mkdir(parents=True, exist_ok=True)
    name = "cnmf"
    a = ad.read_h5ad(src)
    counts = ad.AnnData(X=a.layers["counts"], obs=a.obs[[]].copy(), var=a.var[[]].copy())
    keep = np.asarray(counts.X.sum(axis=0)).ravel() > 0
    counts = counts[:, keep].copy()
    counts_fn = work / "counts.h5ad"
    genes_fn = work / "hvg.txt"
    hvg = [g for g in a.var_names[a.var["highly_variable"].to_numpy(bool)] if g in counts.var_names]
    grid = [int(x) for x in os.environ.get("CNMF_K_GRID", "").split()] or K_GRID
    ks = [k for k in grid if k < len(hvg)]
    n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    if not args.skip_factorize and (work / name).is_dir():
        shutil.rmtree(work / name)    # replicates of an earlier run (other settings) must not be mixed in
    obj = cNMF(output_dir=str(work), name=name)       # (creates its folders)
    if not args.skip_factorize:
        counts.write_h5ad(counts_fn)
        genes_fn.write_text("\n".join(hvg) + "\n")
        log(f"  {counts.n_obs:,} cells x {counts.n_vars:,} genes ({len(hvg)} highly variable); K {ks}, "
            f"{args.n_iter} replicates each, {n_workers} workers")
        obj.prepare(counts_fn=str(counts_fn), components=ks, n_iter=args.n_iter, seed=config.RANDOM_SEED,
                    genes_file=str(genes_fn))
        # one BLAS thread per worker: spawned workers read these when they import numpy
        saved = {v: os.environ.get(v) for v in THREAD_VARS}
        os.environ.update({v: "1" for v in THREAD_VARS})
        try:
            with get_context("spawn").Pool(n_workers) as pool:
                pool.map(_worker, [(str(work), name, i, n_workers) for i in range(n_workers)])
        finally:
            for v, val in saved.items():
                if val is None:
                    os.environ.pop(v, None)
                else:
                    os.environ[v] = val
        obj.combine()
        obj.k_selection_plot(close_fig=True)
    from cnmf.cnmf import load_df_from_npz
    stats = load_df_from_npz(obj.paths["k_selection_stats"])
    k = args.k or choose_k(stats)
    stats["chosen"] = stats.k == k
    log("  K selection: " + ", ".join(f"K={int(r.k)} stability {r.silhouette:.3f}" for r in stats.itertuples())
        + f" -> K={k}")
    man.write(stats, "k_selection", "Per K: stability (silhouette of replicate programmes), prediction error; "
              "chosen K", subdir=SUBDIR)
    obj.consensus(k, density_threshold=DENSITY_THRESHOLD, show_clustering=False, close_clustergram_fig=True)
    usage, scores, _, _ = obj.load_results(K=k, density_threshold=DENSITY_THRESHOLD, n_top_genes=N_TOP)
    factors = [f"P{i + 1:02d}" for i in range(usage.shape[1])]
    usage.columns, scores.columns = factors, factors
    usage = usage.loc[a.obs_names]
    terms = go_terms()
    info = []
    for f in factors:
        topg = list(scores[f].sort_values(ascending=False).index[:N_TOP])
        tg = {g.upper() for g in topg}
        best = max(terms.items(), key=lambda kv: len(tg & kv[1]) / min(len(tg), len(kv[1])),
                   default=("", set()))
        info.append({"factor": f, "scope": "global", "best_prior_set": "", "overlap_coefficient": 0.0,
                     "label": "new", "free": True, "nearest_go": best[0],
                     "nearest_go_overlap": len(tg & best[1]) / min(len(tg), max(len(best[1]), 1)),
                     "top_genes": "|".join(topg[:20]), "mean_cell_score": float(usage[f].mean())})
    info = pd.DataFrame(info)
    man.write(info, "factor_info", "Per programme: top genes, nearest GO BP term, mean usage", subdir=SUBDIR)
    w = scores.T
    w.index.name = "factor"
    man.write(w.reset_index(), "factor_gene_weights", "Programmes x genes: cNMF gene-spectra scores", subdir=SUBDIR)
    write_group_scores(man, a.obs, usage.to_numpy(), factors, SUBDIR)
    man.write(pd.DataFrame([{"k": k, "replicates_per_k": args.n_iter, "k_grid": " ".join(map(str, ks)),
                             "density_threshold": DENSITY_THRESHOLD, "n_cells": counts.n_obs,
                             "n_genes_hvg": len(hvg)}]),
              "fit_summary", "cNMF: chosen K, replicates, cells, genes", subdir=SUBDIR)
    man.flush()


def main() -> None:
    p = cli.build_parser(__doc__)
    p.add_argument("--k", type=int, default=int(os.environ.get("CNMF_K", 0)) or None,
                   help="number of programmes (default: chosen from the stability statistics)")
    p.add_argument("--n-iter", type=int, default=int(os.environ.get("CNMF_N_ITER", N_ITER)),
                   help="NMF replicates per K")
    p.add_argument("--skip-factorize", action="store_true",
                   help="reuse the replicates of an earlier run; only redo the consensus (e.g. for another K)")
    args = p.parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
