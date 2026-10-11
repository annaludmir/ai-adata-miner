#!/usr/bin/env python3
"""27 - Fit scHPF: gene programmes with no prior knowledge (a second method), per stratum.

Single-cell hierarchical Poisson factorisation (Levitin et al., Mol Syst Biol
2019; package scHPF, github.com/simslab/scHPF) models raw counts directly as
Poisson draws from a mix of programmes, with priors that absorb each cell's and
each gene's overall depth. Like cNMF (script 26) no gene set steers it; a second
prior-free method shows which programmes do not depend on the factorisation
chosen.

Runs on the same cells and genes as Spectra and cNMF: the stratified subsample
of script 24 (<AIM_WORK>/spectra/<ns>/input.h5ad, raw counts) and its highly
variable genes. K = the dataset's cNMF K (26_cnmf/fit_summary.csv, after
--finalise the same for both chemistries) so the two prior-free fits have the
same resolution; SCHPF_K overrides it, and without a cNMF fit K = DEFAULT_K.
N_TRIALS random starts (fixed seed), the one with the lowest loss kept.

Outputs
  <AIM_WORK>/schpf/<ns>/model.joblib   the fitted model
  csv_exports/<ns>/27_schpf/
    factor_info.csv             per programme: top genes, nearest GO BP term, mean cell score
                                (scope 'global', label 'new', free: the format of 25_spectra)
    factor_gene_weights.csv     programmes x genes: scHPF gene scores
    factor_scores_by_group.csv  mean cell score per programme x (class | class x age | class x region x age | donor)
    fit_summary.csv             K and where it came from, trials, final loss, cells, genes
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import Manifest, log
from lib.programs import write_group_scores

SCRIPT = "27_schpf"
SUBDIR = "27_schpf"
N_TRIALS = 5
DEFAULT_K = 30
N_TOP = 50


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


def pick_k(ns: str, override: int | None) -> tuple[int, str]:
    if override:
        return override, "SCHPF_K"
    f = config.CSV_EXPORTS / ns / "26_cnmf" / "fit_summary.csv"
    if f.exists():
        return int(pd.read_csv(f).k.iloc[0]), "cNMF K (26_cnmf)"
    return DEFAULT_K, "default (no cNMF fit)"


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    import anndata as ad
    import joblib
    import scipy.sparse as sp
    import schpf
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    man = Manifest(ns, SCRIPT)
    src = config.WORK_DIR / "spectra" / ns / "input.h5ad"
    if not src.exists():
        log(f"  {src} missing -- run script 24 first; skipping")
        man.flush()
        return
    work = config.WORK_DIR / "schpf" / ns
    work.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(src)
    hvg = a.var_names[a.var["highly_variable"].to_numpy(bool)]
    X = sp.csr_matrix(a[:, hvg].layers["counts"]).astype(np.float64)
    expressed = np.asarray(X.sum(axis=0)).ravel() > 0
    keep_cells = np.asarray(X[:, expressed].sum(axis=1)).ravel() > 0
    X, genes, cells = X[keep_cells][:, expressed], hvg[expressed], a.obs_names[keep_cells]
    k, k_from = pick_k(ns, args.k)
    log(f"  {X.shape[0]:,} cells x {X.shape[1]:,} genes; K={k} ({k_from}); {args.ntrials} trials")
    np.random.seed(config.RANDOM_SEED)
    model = schpf.run_trials(X.tocoo(), nfactors=k, ntrials=args.ntrials, verbose=False)
    loss = float(model.loss[-1]) if len(model.loss) else np.nan
    log(f"  best of {args.ntrials} trials: {len(model.loss)} loss checks, final loss {loss:.5g}")
    try:
        joblib.dump(model, work / "model.joblib")
    except Exception as err:          # the CSVs hold the results
        log(f"  model not saved: {err}")
    factors = [f"H{i + 1:02d}" for i in range(k)]
    gscore = pd.DataFrame(model.gene_score(), index=genes, columns=factors)
    cscore = pd.DataFrame(model.cell_score(), index=cells, columns=factors)
    terms = go_terms()
    info = []
    for f in factors:
        topg = list(gscore[f].sort_values(ascending=False).index[:N_TOP])
        tg = {g.upper() for g in topg}
        best = max(terms.items(), key=lambda kv: len(tg & kv[1]) / min(len(tg), len(kv[1])),
                   default=("", set()))
        info.append({"factor": f, "scope": "global", "best_prior_set": "", "overlap_coefficient": 0.0,
                     "label": "new", "free": True, "nearest_go": best[0],
                     "nearest_go_overlap": len(tg & best[1]) / min(len(tg), max(len(best[1]), 1)),
                     "top_genes": "|".join(topg[:20]), "mean_cell_score": float(cscore[f].mean())})
    man.write(pd.DataFrame(info), "factor_info", "Per programme: top genes, nearest GO BP term, mean cell score",
              subdir=SUBDIR)
    w = gscore.T
    w.index.name = "factor"
    man.write(w.reset_index(), "factor_gene_weights", "Programmes x genes: scHPF gene scores", subdir=SUBDIR)
    write_group_scores(man, a.obs.loc[cells], cscore.to_numpy(), factors, SUBDIR)
    man.write(pd.DataFrame([{"k": k, "k_from": k_from, "trials": args.ntrials, "final_loss": loss,
                             "loss_checks": len(model.loss), "n_cells": X.shape[0], "n_genes_hvg": X.shape[1]}]),
              "fit_summary", "scHPF: K and its source, trials, final loss, cells, genes", subdir=SUBDIR)
    man.flush()


def main() -> None:
    p = cli.build_parser(__doc__)
    p.add_argument("--k", type=int, default=int(os.environ.get("SCHPF_K", 0)) or None,
                   help="number of programmes (default: the dataset's cNMF K)")
    p.add_argument("--ntrials", type=int, default=int(os.environ.get("SCHPF_TRIALS", N_TRIALS)),
                   help="random starts; the lowest-loss one is kept")
    args = p.parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
