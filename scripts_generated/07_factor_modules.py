#!/usr/bin/env python3
"""07 - Gene modules from the stored factorisation (varm/Loadings).

human_dev ships a precomputed factorisation: varm['Loadings'] is genes x K and
obsm['Factors'] is cells x K.  That is a ready-made set of gene modules that
cost nothing to extract -- exactly what this repo is named after.

For every factor we export its top-loading genes (the module), measure how
concentrated the module is, check how much modules overlap each other, and test
each module for enrichment in the NDD panels by hypergeometric test.  The last
one is the headline result: which transcriptional programme of the developing
brain is most loaded with autism / ID / epilepsy risk genes.

Outputs (csv_exports/<dataset>/07_modules/)
  module_genes_long.csv        factor -> top genes with loadings (tidy)
  module_summary.csv           per factor: concentration, loading stats
  module_overlap_jaccard.csv   Jaccard between top-gene sets of each factor pair
  module_panel_enrichment.csv  hypergeometric enrichment of each panel per module
  gene_factor_loadings_panel.csv  full loading matrix restricted to panel genes
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import (Manifest, gene_frame, list_h5ad_keys, log,
                          read_elem_at, read_var)
from lib.panels import panel_long_frame
from lib.stats_utils import entropy_specificity, gini

SCRIPT = "07_factor_modules"
SUBDIR = "07_modules"


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    keys = list_h5ad_keys(path)
    if "Loadings" not in keys["varm"]:
        log(f"  no varm/Loadings in this file (varm keys: {keys['varm'] or 'none'}); "
            "gene modules must be derived from pseudobulk instead -- see script 10")
        man.flush()
        return

    loadings = np.asarray(read_elem_at(path, "varm/Loadings"), dtype=np.float64)
    var = read_var(path)
    genes = gene_frame(var, key)
    if loadings.shape[0] != len(genes):
        log(f"  Loadings has {loadings.shape[0]} rows but var has {len(genes)}; aborting")
        man.flush()
        return
    n_genes, K = loadings.shape
    log(f"  Loadings: {n_genes:,} genes x {K} factors")

    symbols = genes["symbol"].astype(str).to_numpy()
    factor_names = [f"Factor{i+1}" for i in range(K)]
    signed = bool((loadings < 0).any())
    log(f"  loadings are {'signed' if signed else 'non-negative'}")

    # -- top genes per factor = the module ----------------------------------
    top_n = min(config.TOP_GENES_PER_FACTOR, n_genes)
    rows = []
    module_sets: dict[str, set[str]] = {}
    for j, fname in enumerate(factor_names):
        col = loadings[:, j]
        order = np.argsort(-col)[:top_n]
        module_sets[fname] = set(symbols[order])
        rows.append(pd.DataFrame({
            "factor": fname,
            "rank": np.arange(1, len(order) + 1),
            "gene": symbols[order],
            "var_index": genes.index.to_numpy()[order],
            "loading": col[order],
            "direction": "positive",
        }))
        if signed:
            order_neg = np.argsort(col)[:top_n]
            rows.append(pd.DataFrame({
                "factor": fname,
                "rank": np.arange(1, len(order_neg) + 1),
                "gene": symbols[order_neg],
                "var_index": genes.index.to_numpy()[order_neg],
                "loading": col[order_neg],
                "direction": "negative",
            }))
    long = pd.concat(rows, ignore_index=True)
    man.write(long, "module_genes_long",
              f"Top {top_n} loading genes per factor -- the gene modules", subdir=SUBDIR)

    # -- how concentrated is each module? -----------------------------------
    absl = np.abs(loadings)
    summary = pd.DataFrame({
        "factor": factor_names,
        "loading_max": loadings.max(axis=0),
        "loading_min": loadings.min(axis=0),
        "loading_mean": loadings.mean(axis=0),
        "loading_sd": loadings.std(axis=0),
        "gini_abs_loading": [gini(absl[:, j]) for j in range(K)],
        "concentration": entropy_specificity(absl.T),
        "n_genes_above_1sd": [(loadings[:, j] > loadings[:, j].mean()
                               + loadings[:, j].std()).sum() for j in range(K)],
        "top_gene": [symbols[np.argmax(loadings[:, j])] for j in range(K)],
    })
    man.write(summary, "module_summary",
              "Per factor: loading distribution and how concentrated the module is",
              subdir=SUBDIR)

    # -- do modules describe distinct gene sets? ----------------------------
    jac = []
    for a, b in combinations(factor_names, 2):
        inter = len(module_sets[a] & module_sets[b])
        union = len(module_sets[a] | module_sets[b])
        jac.append({"factor_a": a, "factor_b": b, "n_shared_genes": inter,
                    "jaccard": inter / union if union else np.nan})
    if jac:
        man.write(pd.DataFrame(jac).sort_values("jaccard", ascending=False),
                  "module_overlap_jaccard",
                  "Jaccard overlap of top-gene sets between factor pairs -- high "
                  "values mean the factors are not independent programmes",
                  subdir=SUBDIR)

    # -- are NDD risk genes concentrated in particular modules? -------------
    from scipy.stats import hypergeom
    panels = panel_long_frame()
    universe = set(symbols)
    enr = []
    for (group, panel), sub in panels.groupby(["panel_group", "panel"]):
        panel_genes = set(sub["gene"]) & universe
        if len(panel_genes) < 3:
            continue
        for fname, mod in module_sets.items():
            hits = mod & panel_genes
            k_obs, M, nn, N = len(hits), len(universe), len(panel_genes), len(mod)
            expected = nn * N / M
            enr.append({
                "panel_group": group, "panel": panel, "factor": fname,
                "n_panel_genes_in_data": nn, "module_size": N,
                "n_overlap": k_obs, "expected_overlap": expected,
                "fold_enrichment": k_obs / expected if expected > 0 else np.nan,
                # survival function at k-1 = P(X >= k)
                "hypergeom_p": float(hypergeom.sf(k_obs - 1, M, nn, N)) if k_obs > 0 else 1.0,
                "overlapping_genes": "|".join(sorted(hits)),
            })
    if enr:
        edf = pd.DataFrame(enr)
        # Benjamini-Hochberg across all module x panel tests
        edf = edf.sort_values("hypergeom_p").reset_index(drop=True)
        m = len(edf)
        edf["fdr_bh"] = np.minimum.accumulate(
            (edf["hypergeom_p"] * m / (edf.index + 1))[::-1])[::-1].clip(0, 1)
        man.write(edf.sort_values(["panel_group", "panel", "hypergeom_p"]),
                  "module_panel_enrichment",
                  "Hypergeometric enrichment of each gene panel in each module, "
                  "with BH-FDR -- which programme carries the NDD risk genes",
                  subdir=SUBDIR)
        sig = edf[(edf["fdr_bh"] < 0.05) & (edf["fold_enrichment"] > 1)]
        log(f"  {len(sig)} module x panel enrichments at FDR < 0.05")

    # -- full loading matrix for panel genes only ---------------------------
    panel_symbols = set(panels["gene"])
    mask = np.array([s in panel_symbols for s in symbols])
    if mask.any():
        pm = pd.DataFrame(loadings[mask], columns=factor_names)
        pm.insert(0, "gene", symbols[mask])
        man.write(pm, "gene_factor_loadings_panel",
                  "Full factor loading matrix restricted to panel genes",
                  subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    if getattr(args, "chemistry", "all") != "all":
        log("  (this step describes the file as a whole; chemistry split not applied)")
    for key in cli.selected_datasets(args):
        run(key, args, None, key)


if __name__ == "__main__":
    main()
