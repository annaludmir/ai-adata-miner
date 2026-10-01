#!/usr/bin/env python3
"""15 - Global relationships between groups: correlation, clustering, PCA.

The step that answers "does the data separate the way the biology says it
should?" before any differential test is run. If samples cluster by donor or by
chemistry rather than by cell class, every downstream result is suspect, and
this is where that shows up.

Hierarchical clustering uses 1 - Pearson distance with all three linkages:
single (chains, permissive), complete (strictest), average (the usual default).
Exporting all three makes a structure that only appears under one of them
visible as the artefact it probably is.

Outputs (csv_exports/<dataset>/15_relationships/)
  correlation_<grouping>_<method>.csv  group x group Pearson / Spearman matrix
  linkage_<grouping>_<linkage>.csv     merge order and heights (dendrogram)
  clusters_<grouping>.csv              membership at several cut levels
  pca_scores_<grouping>.csv            group coordinates on PC1..PCn
  pca_variance_<grouping>.csv          variance explained and cumulative
  pca_loadings_<grouping>.csv          gene contribution to each PC
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.bulk_stats import cut_tree_clusters, hierarchical_linkage, pca
from lib.io_utils import Manifest, log

SCRIPT = "15_sample_relationships"
SUBDIR = "15_relationships"
LINKAGES = ["average", "complete", "single"]
TOP_VARIABLE_GENES = 2000
N_COMPONENTS = 10
MAX_GROUPS_FOR_MATRIX = 400


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.set_index(df.columns[0])
    df.index.name = "gene"
    return df.apply(pd.to_numeric, errors="coerce")


def run(key: str, args) -> None:
    cli.banner(SCRIPT, key)
    man = Manifest(key, SCRIPT)
    pb = config.CSV_EXPORTS / key / "09_pseudobulk"
    if not pb.exists():
        log(f"  {pb} not found -- run 09_pseudobulk.py first")
        man.flush()
        return

    for mpath in sorted(pb.glob("*__mean_lognorm.csv")):
        grouping = mpath.name.replace("__mean_lognorm.csv", "")
        expr = load_matrix(mpath).dropna(how="all")
        if expr.shape[1] < 3:
            log(f"  skipping {grouping}: only {expr.shape[1]} group(s)")
            continue
        if expr.shape[1] > MAX_GROUPS_FOR_MATRIX:
            log(f"  skipping {grouping}: {expr.shape[1]} groups exceeds matrix cap")
            continue

        # Restrict to the most variable genes: invariant genes add no signal to
        # a correlation or a PC, only noise and runtime.
        var = expr.var(axis=1).sort_values(ascending=False)
        sub = expr.loc[var[var > 0].head(TOP_VARIABLE_GENES).index]
        if sub.shape[0] < 10:
            log(f"  skipping {grouping}: only {sub.shape[0]} variable genes")
            continue
        log(f"  {grouping}: {expr.shape[1]} groups, {sub.shape[0]} variable genes")

        for method in ("pearson", "spearman"):
            cm = sub.corr(method=method)
            cm.index.name = "group"
            man.write(cm.reset_index(), f"correlation_{grouping}_{method}",
                      f"{method.capitalize()} correlation between {grouping} groups "
                      f"over the top {sub.shape[0]} variable genes", subdir=SUBDIR)

        cuts = None
        for linkage_method in LINKAGES:
            Z, labels = hierarchical_linkage(sub, method=linkage_method)
            if Z is None:
                continue
            merge = pd.DataFrame(Z, columns=["child_a", "child_b", "distance", "n_in_cluster"])
            merge.insert(0, "merge_step", np.arange(1, len(merge) + 1))
            merge.insert(0, "linkage", linkage_method)
            # Leaves are 0..n-1 and refer to `labels`; internal nodes are n+i.
            merge["child_a_label"] = [labels[int(i)] if int(i) < len(labels) else ""
                                      for i in merge["child_a"]]
            merge["child_b_label"] = [labels[int(i)] if int(i) < len(labels) else ""
                                      for i in merge["child_b"]]
            man.write(merge, f"linkage_{grouping}_{linkage_method}",
                      f"{linkage_method} linkage merge order and heights for "
                      f"{grouping} (dendrogram, leaves 0..n-1 index the labels)",
                      subdir=SUBDIR)
            if linkage_method == "average":
                ks = [k for k in (2, 3, 4, 5, 6, 8, 10) if 2 <= k <= len(labels)]
                cuts = cut_tree_clusters(Z, labels, ks)
        if cuts is not None and not cuts.empty:
            cuts.insert(0, "grouping", grouping)
            man.write(cuts, f"clusters_{grouping}",
                      f"Average-linkage cluster membership of {grouping} groups at "
                      "several cut levels", subdir=SUBDIR)

        scores, loadings, variance = pca(sub, n_components=N_COMPONENTS)
        if scores.empty:
            continue
        variance.insert(0, "grouping", grouping)
        man.write(variance, f"pca_variance_{grouping}",
                  f"Variance explained by each PC of {grouping} (proportion and "
                  "cumulative) -- how many axes the structure really needs",
                  subdir=SUBDIR)
        man.write(scores.reset_index(), f"pca_scores_{grouping}",
                  f"{grouping} group coordinates on PC1..PC{scores.shape[1]}",
                  subdir=SUBDIR)
        top = loadings.reindex(loadings.abs().max(axis=1).sort_values(
            ascending=False).head(1000).index)
        man.write(top.reset_index(), f"pca_loadings_{grouping}",
                  f"Gene loadings on each PC of {grouping} (top 1000 by |loading|) "
                  "-- which genes drive each axis", subdir=SUBDIR)
        log(f"    PC1 {variance['proportion_of_variance'].iloc[0]:.1%}, "
            f"PC1-3 {variance['cumulative_proportion'].iloc[min(2, len(variance)-1)]:.1%}")
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key in cli.selected_datasets(args):
        run(key, args)


if __name__ == "__main__":
    main()
