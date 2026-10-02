#!/usr/bin/env python3
"""16 - Expression patterns: Z-scored gene clustering (K-means + hierarchical).

Reduces thousands of variable genes to a handful of interpretable patterns --
"up in radial glia", "rises with age", "cortex-specific" -- which is what makes
the next functional step tractable at all.

Genes are Z-scored across groups first, and that step is not cosmetic: without
it a highly expressed gene and a lowly expressed one with identical kinetics
cluster apart purely by magnitude, when it is the shape we are after.

K is swept rather than guessed, and both halves of the trade-off are reported:
homogeneity (tightness within a cluster) rises with K while separation
(distinctness between clusters) degrades. The sweep table is how you choose K
honestly instead of defending whichever value you happened to run.

The suggested K is the one with the highest silhouette (correlation distance).
An earlier score, homogeneity minus the max between-centroid r, is kept in the
table but no longer used: the max over K(K-1)/2 centroid pairs rises with K
whatever the data, so it picked K=2 on every real grouping and on synthetic
data with 3, 6 or 9 planted patterns.

Outputs (csv_exports/<dataset>/16_patterns/)
  kmeans_sweep_<grouping>.csv         homogeneity / separation / silhouette vs K
  gene_clusters_<grouping>_k<K>.csv   per-gene cluster membership (K = 6, 9, 12
                                      and the suggested K)
  cluster_profiles_<grouping>_k<K>.csv mean +/- SD Z-score per cluster per group
  cluster_quality_<grouping>_k<K>.csv  per-cluster homogeneity and separation
  zscores_<grouping>.csv               Z-scored matrix (heatmap input)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.bulk_stats import (cluster_quality, kmeans_clusters, silhouette_corr,
                            zscore_rows)
from lib.io_utils import Manifest, log
from lib.panels import panel_long_frame

SCRIPT = "16_expression_patterns"
SUBDIR = "16_patterns"
PRIORITY = ["cell_class", "cell_class_x_age", "cell_class_x_region", "region", "age"]
TOP_VARIABLE_GENES = 3000
K_SWEEP = [2, 3, 4, 6, 8, 9, 12, 16]
K_EXPORT = [6, 9, 12]
MIN_GROUPS = 4
# Silhouette range across K below this means the sweep has no real optimum.
FLAT_SILHOUETTE = 0.05


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
    pb = config.CSV_EXPORTS / ns / "09_pseudobulk"
    if not pb.exists():
        log(f"  {pb} not found -- run 09_pseudobulk.py first")
        man.flush()
        return

    available = {p.name.replace("__mean_lognorm.csv", ""): p
                 for p in sorted(pb.glob("*__mean_lognorm.csv"))}
    chosen = [g for g in PRIORITY if g in available] or list(available)
    panels = panel_long_frame()
    gene_to_panel = panels.groupby("gene")["panel"].apply(
        lambda s: "|".join(sorted(set(s)))).to_dict()

    for grouping in chosen:
        expr = load_matrix(available[grouping]).dropna(how="all")
        if expr.shape[1] < MIN_GROUPS:
            log(f"  skipping {grouping}: {expr.shape[1]} groups (<{MIN_GROUPS})")
            continue
        var = expr.var(axis=1).sort_values(ascending=False)
        sub = expr.loc[var[var > 0].head(TOP_VARIABLE_GENES).index]
        if sub.shape[0] < 50:
            log(f"  skipping {grouping}: only {sub.shape[0]} variable genes")
            continue
        z = zscore_rows(sub)
        log(f"  {grouping}: {z.shape[0]} variable genes x {z.shape[1]} groups")

        man.write(z.round(4).reset_index(), f"zscores_{grouping}",
                  f"Z-scored expression of the top {z.shape[0]} variable genes "
                  f"across {grouping} groups -- the heatmap input", subdir=SUBDIR)

        sweep, fits = [], {}
        for k in [k for k in K_SWEEP if k < z.shape[0]]:
            labels, _, inertia = kmeans_clusters(z, k, seed=config.RANDOM_SEED)
            if len(np.unique(labels)) < 2:
                continue
            q = cluster_quality(z, labels)
            fits[k] = (labels, q)
            homo = float(q["homogeneity_mean_corr_to_centroid"].mean())
            sep = float(q["separation_max_corr_to_other_centroid"].max())
            sweep.append({
                "grouping": grouping, "k": k,
                "mean_homogeneity": homo,
                "min_homogeneity": float(q["homogeneity_mean_corr_to_centroid"].min()),
                "max_separation_corr": sep,
                "mean_separation_corr": float(q["separation_max_corr_to_other_centroid"].mean()),
                # Kept for reference only: biased toward K=2 (see docstring).
                "balance_score": homo - sep,
                "silhouette_corr": silhouette_corr(z, labels),
                "total_within_ss": float(inertia),
                "smallest_cluster": int(q["n_members"].min()),
            })
        if not sweep:
            continue

        sw = pd.DataFrame(sweep)
        best_k = int(sw.loc[sw["silhouette_corr"].idxmax(), "k"])
        sw["suggested"] = sw["k"] == best_k
        man.write(sw, f"kmeans_sweep_{grouping}",
                  "Homogeneity, separation and silhouette across K -- homogeneity "
                  "rises with K while clusters become less distinct; silhouette_corr "
                  "is the suggested criterion, balance_score is biased toward K=2",
                  subdir=SUBDIR)
        best = sw.loc[sw["suggested"]].iloc[0]
        spread = float(sw["silhouette_corr"].max() - sw["silhouette_corr"].min())
        if spread < FLAT_SILHOUETTE:
            log(f"    silhouette is flat across K ({sw['silhouette_corr'].min():.3f}-"
                f"{sw['silhouette_corr'].max():.3f}): no K is clearly preferred; "
                f"K={best_k} exported but read the sweep table, not this pick")
        else:
            log(f"    K with highest silhouette: K={best_k} "
                f"(silhouette {best['silhouette_corr']:.3f}, homogeneity "
                f"{best['mean_homogeneity']:.3f}) "
                f"-- the sweep table is the real answer, this is a starting point")

        # Fixed Ks keep file names stable across runs; the suggested K is
        # added so the pattern the sweep favours is always on disk.
        for k in sorted(set(K_EXPORT) | {best_k}):
            if k not in fits:
                continue
            labels, q = fits[k]
            mem = pd.DataFrame({
                "gene": z.index.astype(str),
                "cluster": labels.astype(int),
            })
            mem["panel"] = mem["gene"].map(gene_to_panel).fillna("")
            mem["peak_group"] = z.columns[np.argmax(z.to_numpy(), axis=1)]
            mem.insert(0, "grouping", grouping)
            man.write(mem.sort_values(["cluster", "gene"]),
                      f"gene_clusters_{grouping}_k{k}",
                      f"Gene -> expression-pattern cluster at K={k} for {grouping}, "
                      "with panel membership and the group where each gene peaks",
                      subdir=SUBDIR)

            # Mean +/- SD profile per cluster: the PDF's cluster plots.
            rows = []
            for c in sorted(np.unique(labels)):
                members = z[labels == c]
                for grp in z.columns:
                    rows.append({
                        "grouping": grouping, "k": k, "cluster": int(c),
                        "n_genes": int(members.shape[0]), "group": grp,
                        "mean_zscore": float(members[grp].mean()),
                        "sd_zscore": float(members[grp].std()),
                        "sem_zscore": float(members[grp].std()
                                            / np.sqrt(max(members.shape[0], 1))),
                    })
            man.write(pd.DataFrame(rows), f"cluster_profiles_{grouping}_k{k}",
                      f"Mean +/- SD Z-score per cluster per group at K={k} -- "
                      "the average expression pattern of each cluster",
                      subdir=SUBDIR)
            q = q.copy()
            q.insert(0, "grouping", grouping)
            q.insert(1, "k", k)
            man.write(q, f"cluster_quality_{grouping}_k{k}",
                      f"Per-cluster homogeneity and separation at K={k}",
                      subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
