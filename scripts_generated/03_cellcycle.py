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
from lib.io_utils import (Manifest, add_derived_obs_columns, log, read_obs,
                          resolve_cluster_columns, resolve_qc_frame, resolve_role)

SCRIPT = "03_cellcycle"
SUBDIR = "03_cellcycle"
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


def run(key: str, args) -> None:
    cli.banner(SCRIPT, key)
    path = cli.resolve_h5ad(key)
    man = Manifest(key, SCRIPT)

    obs = read_obs(path)
    if args.limit_cells:
        obs = obs.iloc[:args.limit_cells]
        log(f"  LIMITED to first {len(obs):,} cells (smoke test)")
    obs = add_derived_obs_columns(obs, key)
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
    for key in cli.selected_datasets(args):
        run(key, args)


if __name__ == "__main__":
    main()
