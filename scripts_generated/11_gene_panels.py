#!/usr/bin/env python3
"""11 - Focused extraction for the NDD and marker gene panels.

The repo's question is where neurodevelopmental-disorder risk genes are
expressed in the developing brain, so the panels get their own tidy export
rather than being buried in a 12,000-row matrix.  Also reports panel coverage,
which matters: a panel gene absent from the file is not an absent gene, and the
distinction has to survive into the downstream analysis.

Marker panels are included so the supplied CellClass labels can be checked
against canonical markers instead of trusted blindly.

Works from the CSVs written by script 09 -- no h5ad access.

Outputs (csv_exports/<dataset>/11_panels/)
  panel_coverage.csv              which panel genes exist in this dataset
  panel_gene_expression_<grp>.csv tidy gene x group expression + detection
  panel_scores_<grouping>.csv     per-group mean panel score, z-scored
  marker_label_check.csv          observed top class for each marker panel
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import Manifest, log
from lib.panels import panel_long_frame

SCRIPT = "11_gene_panels"
SUBDIR = "11_panels"
PRIORITY_GROUPINGS = ["cell_class", "cell_class_x_age", "cell_class_x_region",
                      "region", "age", "sample"]


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.set_index(df.columns[0])
    df.index.name = "gene"
    return df.apply(pd.to_numeric, errors="coerce")


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    man = Manifest(ns, SCRIPT)
    pb_dir = config.CSV_EXPORTS / ns / "09_pseudobulk"
    if not pb_dir.exists():
        log(f"  {pb_dir} not found -- run 09_pseudobulk.py first")
        man.flush()
        return

    panels = panel_long_frame()
    sel_path = pb_dir / "gene_selection.csv"
    if sel_path.exists():
        sel = pd.read_csv(sel_path)
        present = set(sel["symbol"].astype(str))
        exported = set(sel.loc[sel["selected"], "symbol"].astype(str))
    else:
        present = exported = set()

    coverage = panels.copy()
    coverage["present_in_dataset"] = coverage["gene"].isin(present)
    coverage["exported_in_pseudobulk"] = coverage["gene"].isin(exported)
    man.write(coverage, "panel_coverage",
              "Per panel gene: present in this dataset's var, and exported by "
              "script 09 -- absent genes are missing from the file, not silent",
              subdir=SUBDIR)
    summary = coverage.groupby(["panel_group", "panel"]).agg(
        n_genes=("gene", "size"),
        n_present=("present_in_dataset", "sum"),
        n_exported=("exported_in_pseudobulk", "sum")).reset_index()
    summary["fraction_present"] = summary["n_present"] / summary["n_genes"]
    man.write(summary, "panel_coverage_summary",
              "Panel-level coverage counts", subdir=SUBDIR)
    log("  panel coverage: " + ", ".join(
        f"{r.panel}={r.n_present}/{r.n_genes}" for r in summary.itertuples()))

    gene_to_panels = panels.groupby("gene").agg(
        panels=("panel", lambda s: "|".join(sorted(set(s)))),
        panel_groups=("panel_group", lambda s: "|".join(sorted(set(s))))).reset_index()

    matrices = sorted(pb_dir.glob("*__mean_lognorm.csv"))
    groupings = [m.name.replace("__mean_lognorm.csv", "") for m in matrices]
    chosen = [g for g in PRIORITY_GROUPINGS if g in groupings] or groupings
    label_rows = []

    for grouping in chosen:
        expr = load_matrix(pb_dir / f"{grouping}__mean_lognorm.csv")
        dpath = pb_dir / f"{grouping}__detection_fraction.csv"
        detect = load_matrix(dpath) if dpath.exists() else None
        panel_genes = [g for g in expr.index.astype(str) if g in set(panels["gene"])]
        if not panel_genes:
            log(f"  {grouping}: no panel genes in the exported matrix")
            continue
        sub = expr.loc[panel_genes]

        tidy = sub.reset_index().melt(id_vars="gene", var_name="group",
                                      value_name="mean_lognorm")
        if detect is not None:
            dt = detect.reindex(index=panel_genes, columns=expr.columns)
            tidy = tidy.merge(
                dt.reset_index().melt(id_vars="gene", var_name="group",
                                      value_name="detection_fraction"),
                on=["gene", "group"], how="left")
        tidy = tidy.merge(gene_to_panels, on="gene", how="left")
        tidy.insert(0, "grouping", grouping)
        man.write(tidy, f"panel_gene_expression_{grouping}",
                  f"Tidy panel-gene expression and detection per {grouping} group",
                  subdir=SUBDIR)

        # -- per-group panel score, z-scored across groups -------------------
        rows = []
        for (pgroup, pname), g in panels.groupby(["panel_group", "panel"]):
            genes = [x for x in g["gene"] if x in sub.index]
            if len(genes) < 3:
                continue
            mean_expr = sub.loc[genes].mean(axis=0)
            z = ((mean_expr - mean_expr.mean()) / mean_expr.std()
                 if mean_expr.std() > 0 else mean_expr * 0)
            for grp in sub.columns:
                rows.append({"grouping": grouping, "panel_group": pgroup,
                             "panel": pname, "group": grp,
                             "n_genes_used": len(genes),
                             "mean_lognorm": float(mean_expr[grp]),
                             "zscore_across_groups": float(z[grp])})
        if rows:
            scores = pd.DataFrame(rows)
            man.write(scores.sort_values(["panel", "zscore_across_groups"],
                                         ascending=[True, False]),
                      f"panel_scores_{grouping}",
                      f"Mean panel expression per {grouping} group, z-scored across groups",
                      subdir=SUBDIR)
            if grouping == "cell_class":
                for (pg, pn), g in scores.groupby(["panel_group", "panel"]):
                    best = g.loc[g["zscore_across_groups"].idxmax()]
                    label_rows.append({
                        "panel_group": pg, "panel": pn,
                        "highest_scoring_cell_class": best["group"],
                        "zscore": best["zscore_across_groups"],
                        "n_genes_used": int(best["n_genes_used"]),
                    })

    if label_rows:
        check = pd.DataFrame(label_rows)
        # A marker panel should peak on the class it names; anything else is a
        # label-provenance question worth raising before interpreting the data.
        check["expected_match"] = check.apply(
            lambda r: r["panel"].split("_")[0].lower() in
                      str(r["highest_scoring_cell_class"]).lower().replace(" ", "_")
            if r["panel_group"] == "marker" else pd.NA, axis=1)
        man.write(check, "marker_label_check",
                  "Which cell class each panel peaks in -- for marker panels this "
                  "checks the supplied CellClass labels against canonical markers",
                  subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
