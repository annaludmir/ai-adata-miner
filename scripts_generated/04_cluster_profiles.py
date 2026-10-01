#!/usr/bin/env python3
"""04 - Profile every clustering present, and compare them to each other.

Both files ship more than one partition of the same cells (cortex has
Clusters / ClustersModularity / ClustersSurprise / leiden_scVI / louvain).  Each
is profiled against biology (cell class, age, region) *and* against batch
(donor, sample), because a cluster that is pure for one donor is a batch
artefact, not a cell state.  Pairwise ARI/NMI then says whether the alternative
partitions are genuinely different or just relabellings.

Outputs (csv_exports/<dataset>/04_clusters/)
  cluster_profile_<column>.csv   one row per cluster: purity, entropy, batch spread
  cluster_composition_<column>_by_<role>.csv  fractions of each cluster
  clustering_agreement.csv       pairwise ARI / NMI between clusterings
  clustering_overview.csv        one row per clustering column
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
from lib.aggregate import contingency
from lib.io_utils import (Manifest, load_obs, log,
                          resolve_cluster_columns, resolve_qc_frame, resolve_role)
from lib.stats_utils import (adjusted_rand_index, normalized_entropy,
                             normalized_mutual_info)

SCRIPT = "04_cluster_profiles"
SUBDIR = "04_clusters"
PROFILE_ROLES = ["cell_class", "region", "subregion", "cyclephase",
                 "donor", "sample", "chemistry"]
BATCH_ROLES = {"donor", "sample", "chemistry"}
MAX_CLUSTERS = 2000


def purity_and_entropy(obs: pd.DataFrame, cluster_col: str, label_col: str) -> pd.DataFrame:
    """Per cluster: which label dominates, how pure, and how even the mix is."""
    tab = contingency(obs, cluster_col, label_col)
    if tab.empty:
        return pd.DataFrame()
    totals = tab.sum(axis=1)
    dominant = tab.idxmax(axis=1)
    purity = tab.max(axis=1) / totals.replace(0, np.nan)
    ent = tab.apply(lambda r: normalized_entropy(r.to_numpy(dtype=float)), axis=1)
    return pd.DataFrame({
        "cluster": tab.index,
        f"dominant_{label_col}": dominant.values,
        f"purity_{label_col}": purity.values,
        f"entropy_{label_col}": ent.values,
        f"n_distinct_{label_col}": (tab > 0).sum(axis=1).values,
    })


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    obs, _keep = load_obs(path, key, chem, args.limit_cells)
    if obs is None:
        log(f"  no chemistry column in this file -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    if len(obs) == 0:
        log(f"  no cells with chemistry={chem}; skipping")
        man.flush()
        return
    log(f"  {len(obs):,} cells after chemistry filter ({chem or 'pooled'})")
    qc = resolve_qc_frame(obs)
    work = pd.concat([obs, qc], axis=1)
    work = work.loc[:, ~work.columns.duplicated()]

    cluster_cols = resolve_cluster_columns(obs)
    if not cluster_cols:
        log("  no clustering columns found; nothing to profile")
        man.flush()
        return
    log(f"  clusterings: {', '.join(cluster_cols)}")

    roles = {r: resolve_role(obs, r) for r in PROFILE_ROLES}
    roles = {r: c for r, c in roles.items() if c is not None}
    overview = []

    for ccol in cluster_cols:
        n_clusters = obs[ccol].nunique(dropna=True)
        if n_clusters > MAX_CLUSTERS:
            log(f"  skipping {ccol}: {n_clusters:,} clusters exceeds cap")
            continue
        clusters = work[ccol].astype(str)
        base = pd.DataFrame({"cluster": clusters.value_counts().index,
                             "n_cells": clusters.value_counts().values})
        base["fraction_of_dataset"] = base["n_cells"] / len(work)

        profile = base
        for role, col in roles.items():
            pe = purity_and_entropy(work.assign(_c=clusters), "_c", col)
            if pe.empty:
                continue
            pe = pe.rename(columns={
                f"dominant_{col}": f"dominant_{role}",
                f"purity_{col}": f"purity_{role}",
                f"entropy_{col}": f"entropy_{role}",
                f"n_distinct_{col}": f"n_distinct_{role}",
            })
            profile = profile.merge(pe, on="cluster", how="left")

            tab = contingency(work.assign(_c=clusters), "_c", col)
            frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0)
            frac.index.name = "cluster"
            man.write(frac.reset_index(), f"cluster_composition_{ccol}_by_{role}",
                      f"Fraction of each {ccol} cluster made up of each {col}",
                      subdir=SUBDIR)

        if "age_pcw" in work.columns:
            g = work.assign(_c=clusters).groupby("_c", observed=True)["age_pcw"]
            profile = profile.merge(
                pd.DataFrame({"cluster": g.median().index,
                              "age_pcw_median": g.median().values,
                              "age_pcw_min": g.min().values,
                              "age_pcw_max": g.max().values}).assign(
                    age_span_weeks=lambda d: d["age_pcw_max"] - d["age_pcw_min"]),
                on="cluster", how="left")

        for qcc in qc.columns:
            med = work.assign(_c=clusters).groupby("_c", observed=True)[qcc].median()
            profile[f"{qcc}_median"] = profile["cluster"].map(med).astype(float)

        # A cluster that is almost one donor is more likely batch than biology.
        if "purity_donor" in profile.columns:
            profile["batch_dominated_flag"] = profile["purity_donor"] > 0.9
        if "purity_sample" in profile.columns:
            profile["single_sample_flag"] = profile["purity_sample"] > 0.9
        profile["below_min_cells"] = profile["n_cells"] < config.MIN_CELLS_PER_GROUP
        profile.insert(0, "clustering", ccol)
        man.write(profile.sort_values("n_cells", ascending=False),
                  f"cluster_profile_{ccol}",
                  f"Per-cluster profile for {ccol}: purity/entropy vs biology and batch",
                  subdir=SUBDIR)

        row = {"clustering": ccol, "n_clusters": int(n_clusters),
               "median_cluster_size": float(base["n_cells"].median()),
               "min_cluster_size": int(base["n_cells"].min()),
               "max_cluster_size": int(base["n_cells"].max()),
               "n_clusters_below_min_cells": int(profile["below_min_cells"].sum())}
        for flag in ("batch_dominated_flag", "single_sample_flag"):
            if flag in profile.columns:
                row[f"n_{flag}"] = int(profile[flag].sum())
        for role in ("cell_class", "donor"):
            pcol = f"purity_{role}"
            if pcol in profile.columns:
                row[f"mean_{pcol}"] = float(profile[pcol].mean())
        overview.append(row)

    if overview:
        man.write(pd.DataFrame(overview), "clustering_overview",
                  "One row per clustering column: granularity and purity summary",
                  subdir=SUBDIR)

    # -- do the alternative partitions actually disagree? --------------------
    if len(cluster_cols) > 1:
        rows = []
        for a, b in combinations(cluster_cols, 2):
            tab = contingency(obs, a, b)
            if tab.empty:
                continue
            rows.append({
                "clustering_a": a, "clustering_b": b,
                "n_clusters_a": int(obs[a].nunique()), "n_clusters_b": int(obs[b].nunique()),
                "adjusted_rand_index": adjusted_rand_index(tab),
                "normalized_mutual_info": normalized_mutual_info(tab),
            })
        if rows:
            man.write(pd.DataFrame(rows).sort_values("adjusted_rand_index", ascending=False),
                      "clustering_agreement",
                      "Pairwise ARI/NMI between clusterings -- near 1 means the "
                      "partitions are relabellings of the same structure",
                      subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
