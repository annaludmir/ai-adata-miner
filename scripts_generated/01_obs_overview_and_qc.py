#!/usr/bin/env python3
"""01 - Cell-level metadata overview and QC profiling (.obs only).

Answers the questions you must settle before any biology: how many cells per
sample/donor/age, how QC metrics vary across them, and which samples are
technical outliers that could masquerade as biological signal later.

Outputs (csv_exports/<dataset>/01_overview/)
  cell_counts_<col>.csv      cell counts per level of every categorical column
  sample_summary.csv         one row per sample: donor, age, chemistry, QC medians
  donor_summary.csv          one row per donor
  qc_by_<role>.csv           QC distribution per cell class / age / region / phase
  qc_cell_quantiles.csv      global QC distribution (deciles) for threshold setting
  sample_qc_outliers.csv     MAD-based outlier flags per sample per metric
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import (Manifest, add_derived_obs_columns, log, read_obs,
                          resolve_qc_frame, resolve_role)

SCRIPT = "01_obs_overview_and_qc"
SUBDIR = "01_overview"
MAX_LEVELS_FOR_COUNTS = 1000
QUANTILES = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]


def qc_agg(df: pd.DataFrame, by: str, qc_cols: list[str]) -> pd.DataFrame:
    """Per-group QC table: n_cells plus mean/median/sd/IQR of each metric."""
    grouped = df.groupby(by, observed=True)
    out = pd.DataFrame({"n_cells": grouped.size()})
    for col in qc_cols:
        g = grouped[col]
        out[f"{col}_mean"] = g.mean()
        out[f"{col}_median"] = g.median()
        out[f"{col}_sd"] = g.std()
        out[f"{col}_q25"] = g.quantile(0.25)
        out[f"{col}_q75"] = g.quantile(0.75)
    return out.reset_index().rename(columns={by: "group_level"}).assign(grouping=by)


def mad_outlier_flags(df: pd.DataFrame, cols: list[str], n_mads: float = 3.0) -> pd.DataFrame:
    """Robust (median/MAD) z-scores per sample.

    MAD rather than SD because a handful of failed samples would inflate an SD
    enough to hide themselves.
    """
    out = df[["sample"]].copy() if "sample" in df.columns else pd.DataFrame(index=df.index)
    for col in cols:
        v = pd.to_numeric(df[col], errors="coerce")
        med = v.median()
        mad = (v - med).abs().median()
        scale = mad * 1.4826  # MAD -> SD-equivalent for a normal distribution
        z = (v - med) / scale if scale and np.isfinite(scale) and scale > 0 else pd.Series(np.nan, index=v.index)
        out[f"{col}_robust_z"] = z
        out[f"{col}_is_outlier"] = z.abs() > n_mads
    flag_cols = [c for c in out.columns if c.endswith("_is_outlier")]
    out["n_metrics_flagged"] = out[flag_cols].sum(axis=1) if flag_cols else 0
    return out


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
    qc_cols = list(qc.columns)
    log(f"  {len(obs):,} cells | QC metrics available: {', '.join(qc_cols) or 'none'}")
    work = pd.concat([obs, qc], axis=1)
    work = work.loc[:, ~work.columns.duplicated()]

    # -- cell counts per level of every manageable categorical column --------
    for col in obs.columns:
        s = obs[col]
        if pd.api.types.is_numeric_dtype(s) and not isinstance(s.dtype, pd.CategoricalDtype):
            continue
        n_levels = s.nunique(dropna=True)
        if n_levels < 1 or n_levels > MAX_LEVELS_FOR_COUNTS:
            continue
        vc = s.value_counts(dropna=False).rename_axis("level").reset_index(name="n_cells")
        vc["fraction"] = vc["n_cells"] / len(obs)
        vc.insert(0, "column", col)
        man.write(vc, f"cell_counts_{col}".replace("/", "_"),
                  f"Cell counts per level of obs['{col}']", subdir=SUBDIR)

    # -- sample- and donor-level tables -------------------------------------
    sample_col = resolve_role(obs, "sample")
    donor_col = resolve_role(obs, "donor")
    context_roles = ["donor", "age", "chemistry", "region", "subregion",
                     "dissection", "sex", "dev_stage", "transcriptome"]

    for level_name, level_col in (("sample", sample_col), ("donor", donor_col)):
        if level_col is None:
            log(f"  no {level_name} column present; skipping {level_name}_summary")
            continue
        grouped = work.groupby(level_col, observed=True)
        summary = pd.DataFrame({"n_cells": grouped.size()})
        for role in context_roles:
            c = resolve_role(obs, role)
            if c is None or c == level_col:
                continue
            # one representative value + how many distinct ones (>1 = mixed)
            summary[role] = grouped[c].agg(lambda s: str(s.iloc[0]) if len(s) else "")
            summary[f"{role}_n_distinct"] = grouped[c].nunique()
        if "age_pcw" in work.columns:
            summary["age_pcw"] = grouped["age_pcw"].median()
        for qcc in qc_cols:
            summary[f"{qcc}_median"] = grouped[qcc].median()
            summary[f"{qcc}_mean"] = grouped[qcc].mean()
        summary = summary.reset_index().rename(columns={level_col: level_name})
        summary["low_cell_count"] = summary["n_cells"] < config.MIN_CELLS_PER_GROUP
        man.write(summary, f"{level_name}_summary",
                  f"One row per {level_name}: context, cell count and QC medians",
                  subdir=SUBDIR)

        if level_name == "sample" and qc_cols:
            med_cols = [f"{c}_median" for c in qc_cols if f"{c}_median" in summary.columns]
            flags = mad_outlier_flags(summary[["sample"] + med_cols], med_cols)
            flags.insert(1, "n_cells", summary["n_cells"].values)
            man.write(flags, "sample_qc_outliers",
                      "Robust (MAD) z-scores per sample per QC metric; |z|>3 flagged",
                      subdir=SUBDIR)

    # -- QC distribution by biological grouping ------------------------------
    if qc_cols:
        for role in ["cell_class", "region", "subregion", "chemistry", "cyclephase", "donor"]:
            col = resolve_role(obs, role)
            if col is None or obs[col].nunique() > MAX_LEVELS_FOR_COUNTS:
                continue
            man.write(qc_agg(work, col, qc_cols), f"qc_by_{role}",
                      f"QC metric distribution per level of {col}", subdir=SUBDIR)
        if "age_pcw" in work.columns and work["age_pcw"].notna().any():
            man.write(qc_agg(work, "age_pcw", qc_cols), "qc_by_age",
                      "QC metric distribution per post-conception week", subdir=SUBDIR)

        quant = work[qc_cols].quantile(QUANTILES).reset_index().rename(columns={"index": "quantile"})
        man.write(quant, "qc_cell_quantiles",
                  "Global per-cell QC quantiles -- use to set filtering thresholds",
                  subdir=SUBDIR)

    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key in cli.selected_datasets(args):
        run(key, args)


if __name__ == "__main__":
    main()
