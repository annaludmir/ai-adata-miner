#!/usr/bin/env python3
"""18 - What is comparable between the v2 and v3 chemistries?

Runs POOLED by design: it is the one analysis that has to see both chemistries
at once, because its whole purpose is to measure their overlap.

Chemistry is not a nuisance covariate in these datasets, it is a near-complete
confound with age. In human_dev the v2 libraries cover roughly 6-10 pcw and the
v3 libraries 5-5.5 and 11.5-14, so a gene that "rises with development" in the
pooled data may simply be a gene that v3 captures better. Stratifying (the
pipeline default, --chemistry each) removes the confound but costs coverage:
neither chemistry spans development on its own.

This quantifies that trade-off before any of it is interpreted -- which ages,
regions, cell classes and donors exist in both chemistries, and therefore which
contrasts survive stratification and which are simply unavailable.

Outputs (csv_exports/<dataset>/18_chemistry/)
  chemistry_cell_counts.csv      cells, samples and donors per chemistry
  comparability_<role>.csv       per level: n_v2, n_v3, whether both are present
  age_chemistry_overlap.csv      the age axis specifically, with a verdict
  chemistry_summary.csv          one row per dataset: how much is shared
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
                          resolve_role)

SCRIPT = "18_chemistry_comparability"
SUBDIR = "18_chemistry"
ROLES = ["cell_class", "region", "subregion", "dissection", "donor", "sample",
         "cyclephase", "dev_stage"]
MAX_LEVELS = 400
MIN_CELLS = config.MIN_CELLS_PER_GROUP


def comparability(obs: pd.DataFrame, level_col: str, chem_col: str,
                  role: str) -> pd.DataFrame:
    """Per level of `level_col`: how many cells of each chemistry, and a verdict."""
    tab = pd.crosstab(obs[level_col].astype(str), obs[chem_col].astype(str))
    for c in config.CHEMISTRIES:
        if c not in tab.columns:
            tab[c] = 0
    tab = tab[list(config.CHEMISTRIES)]
    out = tab.reset_index()
    out.columns = ["level"] + [f"n_{c}" for c in config.CHEMISTRIES]
    out.insert(0, "role", role)
    n2, n3 = out["n_v2"], out["n_v3"]
    out["total"] = n2 + n3
    out["in_both"] = (n2 > 0) & (n3 > 0)
    # "Both present" is not the same as "both usable": a level with 4 cells in v3
    # cannot support a stratified comparison even though it is non-zero.
    out["comparable"] = (n2 >= MIN_CELLS) & (n3 >= MIN_CELLS)
    out["verdict"] = np.where(
        out["comparable"], "comparable across chemistry",
        np.where(out["in_both"], f"present in both but under {MIN_CELLS} cells in one",
                 np.where(n2 > 0, "v2 only", "v3 only")))
    return out.sort_values("total", ascending=False)


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, None)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    obs = add_derived_obs_columns(read_obs(path), key)
    if args.limit_cells:
        obs = obs.iloc[:args.limit_cells]
    chem_col = resolve_role(obs, "chemistry")
    if chem_col is None:
        log("  no chemistry column in this file; nothing to compare")
        man.flush()
        return
    present = sorted(obs[chem_col].astype(str).unique())
    log(f"  chemistries present: {', '.join(present)}")
    if len(present) < 2:
        log("  only one chemistry in this file -- stratification is a no-op here")

    # -- how much data does each chemistry carry? --------------------------
    rows = []
    sample_col, donor_col = resolve_role(obs, "sample"), resolve_role(obs, "donor")
    for c in present:
        sub = obs[obs[chem_col].astype(str) == c]
        rows.append({
            "dataset": key, "chemistry": c,
            "n_cells": len(sub),
            "fraction_of_cells": len(sub) / len(obs),
            "n_samples": int(sub[sample_col].nunique()) if sample_col else np.nan,
            "n_donors": int(sub[donor_col].nunique()) if donor_col else np.nan,
            "age_min": float(sub["age_pcw"].min()) if "age_pcw" in sub else np.nan,
            "age_max": float(sub["age_pcw"].max()) if "age_pcw" in sub else np.nan,
            "n_ages": int(sub["age_pcw"].nunique()) if "age_pcw" in sub else np.nan,
        })
    counts = pd.DataFrame(rows)
    man.write(counts, "chemistry_cell_counts",
              "Cells, samples, donors and age span carried by each chemistry",
              subdir=SUBDIR)
    for r in counts.itertuples():
        log(f"    {r.chemistry}: {r.n_cells:,} cells ({r.fraction_of_cells:.1%}), "
            f"ages {r.age_min}-{r.age_max}, {r.n_donors} donors")

    # -- the age axis, where the confound actually bites -------------------
    summary = {"dataset": key, "n_cells": len(obs),
               "chemistries": "|".join(present)}
    if "age_pcw" in obs.columns and obs["age_pcw"].notna().any():
        age = comparability(obs.assign(_a=obs["age_pcw"].astype(str)),
                            "_a", chem_col, "age")
        age["age_pcw"] = pd.to_numeric(age["level"], errors="coerce")
        age = age.sort_values("age_pcw")
        man.write(age, "age_chemistry_overlap",
                  "Cells per age per chemistry, with whether that age can support "
                  "a within-age comparison across chemistry", subdir=SUBDIR)
        shared_cells = int(age.loc[age["comparable"], "total"].sum())
        summary.update({
            "n_ages": int(len(age)),
            "n_ages_in_both": int(age["in_both"].sum()),
            "n_ages_comparable": int(age["comparable"].sum()),
            "frac_cells_at_comparable_ages": shared_cells / len(obs),
        })
        log(f"  ages: {len(age)} total, {int(age['in_both'].sum())} in both "
            f"chemistries, {int(age['comparable'].sum())} with >= {MIN_CELLS} "
            f"cells in each ({shared_cells / len(obs):.1%} of cells)")
        if int(age["comparable"].sum()) == 0:
            log("  WARNING: no age is measured by both chemistries. Age and "
                "chemistry are fully confounded -- a developmental trend across "
                "the whole range cannot be separated from the v2->v3 switch.")

    # -- every other covariate ---------------------------------------------
    for role in ROLES:
        col = resolve_role(obs, role)
        if col is None or obs[col].nunique(dropna=True) > MAX_LEVELS:
            continue
        cmp_df = comparability(obs, col, chem_col, role)
        man.write(cmp_df, f"comparability_{role}",
                  f"Per {col} level: cells of each chemistry and whether the "
                  "level survives stratification", subdir=SUBDIR)
        summary[f"{role}_n_levels"] = int(len(cmp_df))
        summary[f"{role}_n_comparable"] = int(cmp_df["comparable"].sum())
        if role == "donor":
            # A donor sequenced on only one chemistry means donor is nested in
            # chemistry, so a donor-level random effect cannot separate them.
            nested = not bool(cmp_df["in_both"].any())
            summary["donor_nested_in_chemistry"] = nested
            if nested:
                log("  every donor is sequenced on a single chemistry: donor is "
                    "nested in chemistry, so the two cannot be separated by "
                    "treating donor as the replicate")

    man.write(pd.DataFrame([summary]), "chemistry_summary",
              "One row per dataset: how much of it survives chemistry stratification",
              subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    # Pooled by construction: comparing the chemistries requires seeing both.
    for key in cli.selected_datasets(args):
        run(key, args, None, key)


if __name__ == "__main__":
    main()
