#!/usr/bin/env python3
"""03 - Cell cycle and proliferation dynamics.

Proliferation is the dominant axis of variation in fetal brain, so it gets its
own extraction: phase composition and continuous cycling scores per cell class,
per week, per region and per cluster.  The age x cell-class table is the
proliferation trajectory -- when each progenitor class stops dividing.

Outputs (csv_exports/<dataset>/03_cellcycle/)
  phase_fractions_by_<grouping>.csv   cycle-phase composition per group
  cycle_scores_by_<grouping>.csv      continuous score stats per group
  proliferation_trajectory.csv        cell class x week: cycling fraction + scores
  phase_fractions_by_cell_class_x_{age,donor}.csv  phase composition per class x week / donor
  phase_fractions_by_cell_class_x_age_x_depth.csv  the same per class x week x UMI quintile
  cycling_summary.csv                 one row per cell class, overall
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.aggregate import contingency
from lib.io_utils import (Manifest, load_obs, log,
                          resolve_cluster_columns, resolve_qc_frame, resolve_role)

SCRIPT = "03_cellcycle"
SUBDIR = "03_cellcycle"
DEPTH_BINS = 5      # UMI quintiles within each cell class
CYCLE_SCORE_COLS = ["cell_cycle_score", "cycling_score", "cc_g1", "cc_s", "cc_g2m"]
MAX_GROUP_LEVELS = 1000


def score_stats(df: pd.DataFrame, by: str, cols: list[str]) -> pd.DataFrame:
    g = df.groupby(by, observed=True)
    out = pd.DataFrame({"n_cells": g.size()})
    for c in cols:
        out[f"{c}_mean"] = g[c].mean()
        out[f"{c}_median"] = g[c].median()
        out[f"{c}_sd"] = g[c].std()
    return out.reset_index().rename(columns={by: "group_level"})


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

    phase_col = "cyclephase_h" if "cyclephase_h" in work.columns else None
    score_cols = [c for c in CYCLE_SCORE_COLS if c in work.columns]
    log(f"  phase column: {phase_col} | score columns: {', '.join(score_cols) or 'none'}")

    # A boolean 'is cycling' from whichever source the file provides.
    if "Cycling" in obs.columns:
        work["is_cycling"] = obs["Cycling"].astype(bool)
    elif "IsCycling" in obs.columns:
        work["is_cycling"] = obs["IsCycling"].astype(float) > 0
    elif phase_col is not None:
        work["is_cycling"] = ~work[phase_col].astype(str).isin(["Non-cycling", "Post-M"])
    else:
        work["is_cycling"] = np.nan

    groupings: dict[str, str] = {}
    for role in ["cell_class", "region", "subregion", "chemistry", "donor", "sample"]:
        col = resolve_role(obs, role)
        if col is not None and obs[col].nunique(dropna=True) <= MAX_GROUP_LEVELS:
            groupings[role] = col
    if "age_pcw" in work.columns:
        groupings["age"] = "age_pcw"
    for col in resolve_cluster_columns(obs):
        if obs[col].nunique(dropna=True) <= MAX_GROUP_LEVELS:
            groupings[f"cluster_{col}"] = col

    for role, col in groupings.items():
        if phase_col is not None:
            counts = contingency(work, col, phase_col)
            if not counts.empty:
                frac = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
                frac.index.name = role
                out = frac.reset_index()
                out.insert(1, "n_cells", counts.sum(axis=1).values)
                man.write(out, f"phase_fractions_by_{role}",
                          f"Cell-cycle phase composition per {col}", subdir=SUBDIR)
        if score_cols:
            stats = score_stats(work, col, score_cols)
            cyc = work.groupby(col, observed=True)["is_cycling"].mean()
            stats["fraction_cycling"] = stats["group_level"].map(cyc).astype(float)
            stats.insert(0, "grouping", role)
            man.write(stats, f"cycle_scores_by_{role}",
                      f"Cell-cycle score statistics per {col}", subdir=SUBDIR)

    # -- the trajectory: cell class x post-conception week -------------------
    class_col = resolve_role(obs, "cell_class")
    if class_col is not None and "age_pcw" in work.columns:
        g = work.groupby([class_col, "age_pcw"], observed=True)
        traj = pd.DataFrame({"n_cells": g.size(), "fraction_cycling": g["is_cycling"].mean()})
        for c in score_cols:
            traj[f"{c}_mean"] = g[c].mean()
            traj[f"{c}_median"] = g[c].median()
        traj = traj.reset_index().rename(columns={class_col: "cell_class"})
        traj["below_min_cells"] = traj["n_cells"] < config.MIN_CELLS_PER_GROUP
        man.write(traj.sort_values(["cell_class", "age_pcw"]), "proliferation_trajectory",
                  "Cycling fraction and cycle scores per cell class per week -- "
                  "shows when each progenitor class exits the cycle", subdir=SUBDIR)

        # Phase composition per cell class x week (and x donor): the input for
        # cell-cycle kinetics within a progenitor type, e.g. whether cycling
        # progenitors spend a growing share of their time in G1 as
        # neurogenesis proceeds.
        if phase_col is not None:
            for by, label in (("age_pcw", "age"), (resolve_role(obs, "donor"), "donor")):
                if by is None or by not in work.columns:
                    continue
                ct = (work.groupby([class_col, by, phase_col], observed=True).size()
                      .unstack(phase_col, fill_value=0))
                fr = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)
                fr.insert(0, "n_cells", ct.sum(axis=1))
                fr = fr.reset_index().rename(columns={class_col: "cell_class", by: label})
                fr["below_min_cells"] = fr["n_cells"] < config.MIN_CELLS_PER_GROUP
                man.write(fr, f"phase_fractions_by_cell_class_x_{label}",
                          f"Cell-cycle phase composition per cell class x {label}", subdir=SUBDIR)

            # The same per class x week, within depth bins. Phase calls come from
            # marker expression, so they can shift with UMIs per cell -- which
            # falls with age in cortex. Bins are UMI quintiles within each cell
            # class (all ages of this stratum), so a downstream analysis can
            # compare ages at matched depth.
            if "total_umis" in work.columns and work["total_umis"].notna().any():
                d = work[[class_col, "age_pcw", phase_col, "total_umis"]].dropna()
                d = d.assign(depth_bin=d.groupby(class_col, observed=True)["total_umis"].transform(
                    lambda x: pd.qcut(x.rank(method="first"), DEPTH_BINS, labels=False)
                    if len(x) >= DEPTH_BINS else 0))
                ct = (d.groupby([class_col, "age_pcw", "depth_bin", phase_col], observed=True).size()
                      .unstack(phase_col, fill_value=0))
                fr = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)
                fr.insert(0, "n_cells", ct.sum(axis=1))
                med = d.groupby([class_col, "age_pcw", "depth_bin"], observed=True)["total_umis"].median()
                fr.insert(1, "median_umis", med.reindex(fr.index).to_numpy())
                fr = fr.reset_index().rename(columns={class_col: "cell_class", "age_pcw": "age"})
                man.write(fr, "phase_fractions_by_cell_class_x_age_x_depth",
                          "Cell-cycle phase composition per cell class x age x UMI quintile "
                          "(quintiles within each class) -- for depth-matched phase shares",
                          subdir=SUBDIR)

        summary = work.groupby(class_col, observed=True).agg(
            n_cells=("is_cycling", "size"),
            fraction_cycling=("is_cycling", "mean"),
            age_min=("age_pcw", "min"), age_max=("age_pcw", "max"),
            age_median=("age_pcw", "median")).reset_index().rename(
            columns={class_col: "cell_class"})
        man.write(summary, "cycling_summary",
                  "Per cell class: overall cycling fraction and age span", subdir=SUBDIR)

    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
