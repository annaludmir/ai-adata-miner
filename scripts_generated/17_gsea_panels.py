#!/usr/bin/env python3
"""17 - GSEA of the gene panels against each group's ranked gene list.

Complements script 07/11's hypergeometric tests, and exists because of their
known weakness: a hypergeometric test needs a thresholded target set, so a
panel that shifts coherently without any single gene clearing the cutoff is
invisible to it. GSEA walks the *whole* ranking and detects exactly that case.

Each group's genes are ranked by specificity-weighted fold change (from script
10), a panel's members are located in that ranking, and the weighted running
sum gives the enrichment score. Significance comes from permuting gene sets of
the same size.

Caveat carried into the output: gene-set permutation does not preserve
gene-gene correlation, so p-values are anti-conservative for co-regulated
panels. The ranking of hypotheses is reliable; the absolute p is optimistic.

Outputs (csv_exports/<dataset>/17_gsea/)
  gsea_results.csv        ES / NES / p / BH-FDR per (group, panel)
  gsea_leading_edge.csv   the core genes driving each significant enrichment
  gsea_running_sum_<...>  running-sum curve for the top enrichments
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.bulk_stats import benjamini_hochberg, gsea_enrichment_score, gsea_test
from lib.io_utils import Manifest, log
from lib.panels import panel_long_frame

SCRIPT = "17_gsea_panels"
SUBDIR = "17_gsea"
PRIORITY = ["cell_class", "cell_class_x_age", "region", "age"]
N_PERMUTATIONS = 1000
MIN_SET_SIZE = 5
MAX_GROUPS = 60
TOP_CURVES = 10


def run(key: str, args) -> None:
    cli.banner(SCRIPT, key)
    man = Manifest(key, SCRIPT)
    mk = config.CSV_EXPORTS / key / "10_markers"
    if not mk.exists():
        log(f"  {mk} not found -- run 10_marker_specificity.py first")
        man.flush()
        return

    panels = panel_long_frame()
    results, leading_rows, curves = [], [], []

    available = {p.name.replace("specificity_", "").replace(".csv", ""): p
                 for p in sorted(mk.glob("specificity_*.csv"))}
    chosen = [g for g in PRIORITY if g in available] or list(available)[:2]

    for grouping in chosen:
        spec = pd.read_csv(available[grouping])
        needed = {"gene", "top_group", "log2fc_vs_rest"}
        if not needed <= set(spec.columns):
            log(f"  skipping {grouping}: missing {needed - set(spec.columns)}")
            continue
        groups = spec["top_group"].dropna().unique()
        if len(groups) > MAX_GROUPS:
            log(f"  skipping {grouping}: {len(groups)} groups exceeds cap {MAX_GROUPS}")
            continue
        log(f"  {grouping}: {len(groups)} groups x {spec['gene'].nunique():,} genes")

        # Rank every gene once per group: genes assigned to this group keep
        # their fold change, the rest are ranked below by construction.
        for grp in groups:
            score = np.where(spec["top_group"] == grp,
                             spec["log2fc_vs_rest"].fillna(0.0),
                             -np.abs(spec["log2fc_vs_rest"].fillna(0.0)))
            if "tau_specificity" in spec.columns:
                score = score * spec["tau_specificity"].fillna(0.0).to_numpy()
            order = np.argsort(-score)
            ranked = spec["gene"].astype(str).to_numpy()[order].tolist()
            weights = np.abs(score[order])

            for (pgroup, pname), g in panels.groupby(["panel_group", "panel"]):
                members = set(g["gene"]) & set(ranked)
                if len(members) < MIN_SET_SIZE:
                    continue
                res = gsea_test(ranked, members, weights=weights,
                                n_permutations=N_PERMUTATIONS,
                                seed=config.RANDOM_SEED)
                if not np.isfinite(res["es"]):
                    continue
                results.append({
                    "grouping": grouping, "group": str(grp),
                    "panel_group": pgroup, "panel": pname,
                    "n_genes_in_set": res["n_genes_in_set"],
                    "es": res["es"], "nes": res["nes"],
                    "p_value": res["p_value"],
                    "peak_rank": res["peak_rank"],
                    "leading_edge_size": res.get("leading_edge_size", 0),
                    "leading_edge_genes": "|".join(res["leading_edge"][:60]),
                })
                for gene in res["leading_edge"]:
                    leading_rows.append({
                        "grouping": grouping, "group": str(grp),
                        "panel_group": pgroup, "panel": pname, "gene": gene,
                        "es": res["es"], "p_value": res["p_value"],
                    })

    if not results:
        log("  no panel reached the minimum set size in any ranking")
        man.flush()
        return

    df = pd.DataFrame(results)
    df["fdr_bh"] = benjamini_hochberg(df["p_value"].to_numpy())
    df["significant"] = (df["fdr_bh"] < 0.05) & (df["es"] > 0)
    df = df.sort_values(["fdr_bh", "p_value"])
    man.write(df, "gsea_results",
              "GSEA per (group, panel): enrichment score, NES, permutation p and "
              "BH-FDR. Positive ES = panel concentrated among that group's top genes.",
              subdir=SUBDIR)
    log(f"  {int(df['significant'].sum())} / {len(df)} enrichments at FDR < 0.05")

    if leading_rows:
        led = pd.DataFrame(leading_rows)
        man.write(led, "gsea_leading_edge",
                  "Leading-edge genes: the core subset driving each enrichment "
                  "signal -- the genes to actually follow up", subdir=SUBDIR)

    # Running-sum curves for the strongest results, so the walk is inspectable.
    top = df[df["significant"]].head(TOP_CURVES)
    if top.empty:
        top = df.head(min(3, len(df)))
    for _, row in top.iterrows():
        sp = pd.read_csv(mk / f"specificity_{row['grouping']}.csv")
        score = np.where(sp["top_group"].astype(str) == row["group"],
                         sp["log2fc_vs_rest"].fillna(0.0),
                         -np.abs(sp["log2fc_vs_rest"].fillna(0.0)))
        if "tau_specificity" in sp.columns:
            score = score * sp["tau_specificity"].fillna(0.0).to_numpy()
        order = np.argsort(-score)
        ranked = sp["gene"].astype(str).to_numpy()[order].tolist()
        members = set(panels.loc[panels["panel"] == row["panel"], "gene"]) & set(ranked)
        es, running, _ = gsea_enrichment_score(ranked, members, np.abs(score[order]))
        if not np.isfinite(es):
            continue
        curves.append(pd.DataFrame({
            "grouping": row["grouping"], "group": row["group"],
            "panel": row["panel"], "rank": np.arange(1, len(running) + 1),
            "running_sum": running,
            "in_set": [g in members for g in ranked],
        }))
    if curves:
        man.write(pd.concat(curves, ignore_index=True), "gsea_running_sum_top",
                  "Running-sum curves for the strongest enrichments; ES is the "
                  "maximum deviation from zero along the walk", subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key in cli.selected_datasets(args):
        run(key, args)


if __name__ == "__main__":
    main()
