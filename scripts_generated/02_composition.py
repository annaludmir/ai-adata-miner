#!/usr/bin/env python3
"""02 - Cellular composition across every biological and technical grouping.

For each (grouping x label) pair this writes the raw contingency table, the
row-normalised fractions, and two enrichment views -- log2(observed/expected)
and standardised chi-square residuals -- so a composition shift can be read off
directly instead of eyeballed from counts.  Per-group diversity indices and a
Cramer's V summary say which groupings actually restructure composition.

Outputs (csv_exports/<dataset>/02_composition/)
  counts_<label>_by_<grouping>.csv      cells per (group, label)
  fractions_<label>_by_<grouping>.csv   row-normalised composition
  log2oe_<label>_by_<grouping>.csv      log2(observed/expected) enrichment
  residuals_<label>_by_<grouping>.csv   standardised chi-square residuals
  diversity_<label>_by_<grouping>.csv   entropy / Simpson / dominant label
  composition_long.csv                  tidy long form of every pair
  association_summary.csv               Cramer's V per (grouping, label) pair
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
                          resolve_cluster_columns, resolve_role)
from lib.stats_utils import (chi2_standardised_residuals, cramers_v,
                             log2_observed_expected, normalized_entropy,
                             shannon_entropy, simpson_index)

SCRIPT = "02_composition"
SUBDIR = "02_composition"
GROUPING_ROLES = ["age", "region", "subregion", "dissection", "donor", "sample",
                  "chemistry", "dev_stage", "age_chem"]
LABEL_ROLES = ["cell_class", "cyclephase"]
MAX_GROUP_LEVELS = 1000


def diversity_table(counts: pd.DataFrame) -> pd.DataFrame:
    """Per-row composition diversity plus which label dominates."""
    rows = []
    for level, row in counts.iterrows():
        v = row.to_numpy(dtype=float)
        total = v.sum()
        rows.append({
            "group_level": level,
            "n_cells": int(total),
            "n_labels_observed": int((v > 0).sum()),
            "shannon_entropy": shannon_entropy(v),
            "normalized_entropy": normalized_entropy(v),
            "simpson_index": simpson_index(v),
            "dominant_label": counts.columns[int(np.argmax(v))] if total > 0 else "",
            "dominant_fraction": float(v.max() / total) if total > 0 else np.nan,
            "below_min_cells": bool(total < config.MIN_CELLS_PER_GROUP),
        })
    return pd.DataFrame(rows)


def run(key: str, args) -> None:
    cli.banner(SCRIPT, key)
    path = cli.resolve_h5ad(key)
    man = Manifest(key, SCRIPT)

    obs = read_obs(path)
    if args.limit_cells:
        obs = obs.iloc[:args.limit_cells]
        log(f"  LIMITED to first {len(obs):,} cells (smoke test)")
    obs = add_derived_obs_columns(obs, key)

    # Resolve grouping columns; 'age' prefers the derived numeric age_pcw.
    groupings: dict[str, str] = {}
    for role in GROUPING_ROLES:
        col = "age_pcw" if role == "age" and "age_pcw" in obs.columns else resolve_role(obs, role)
        if col is None or col not in obs.columns:
            continue
        if obs[col].nunique(dropna=True) > MAX_GROUP_LEVELS:
            log(f"  skipping grouping '{role}' ({col}): {obs[col].nunique():,} levels")
            continue
        groupings[role] = col
    for col in resolve_cluster_columns(obs):
        if obs[col].nunique(dropna=True) <= MAX_GROUP_LEVELS:
            groupings[f"cluster_{col}"] = col

    labels: dict[str, str] = {}
    for role in LABEL_ROLES:
        col = "cyclephase_h" if role == "cyclephase" and "cyclephase_h" in obs.columns \
            else resolve_role(obs, role)
        if col is not None and col in obs.columns:
            labels[role] = col

    log(f"  groupings: {', '.join(groupings)}")
    log(f"  labels:    {', '.join(labels)}")
    if not labels:
        log("  no label column (cell class / cycle phase) found -- nothing to do")
        man.flush()
        return

    long_rows, assoc_rows = [], []
    for lrole, lcol in labels.items():
        for grole, gcol in groupings.items():
            if gcol == lcol:
                continue
            counts = contingency(obs, gcol, lcol)
            if counts.empty or counts.shape[1] < 2:
                continue
            counts.index.name = grole
            tag = f"{lrole}_by_{grole}"

            fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
            man.write(counts.reset_index(), f"counts_{tag}",
                      f"Cell counts: {gcol} x {lcol}", subdir=SUBDIR)
            man.write(fractions.reset_index(), f"fractions_{tag}",
                      f"Row-normalised composition of {lcol} within each {gcol}",
                      subdir=SUBDIR)
            man.write(log2_observed_expected(counts).reset_index(), f"log2oe_{tag}",
                      f"log2(observed/expected) enrichment of {lcol} per {gcol}",
                      subdir=SUBDIR)
            man.write(chi2_standardised_residuals(counts).reset_index(), f"residuals_{tag}",
                      f"Standardised chi-square residuals ({lcol} x {gcol}); |r|>2 notable",
                      subdir=SUBDIR)
            div = diversity_table(counts)
            div.insert(0, "grouping", grole)
            man.write(div, f"diversity_{tag}",
                      f"Composition diversity of {lcol} within each {gcol}", subdir=SUBDIR)

            melted = counts.reset_index().melt(id_vars=grole, var_name="label",
                                               value_name="n_cells")
            melted.columns = ["group_level", "label", "n_cells"]
            melted.insert(0, "grouping", grole)
            melted.insert(1, "label_type", lrole)
            totals = melted.groupby("group_level")["n_cells"].transform("sum")
            melted["fraction"] = melted["n_cells"] / totals.replace(0, np.nan)
            long_rows.append(melted)

            assoc_rows.append({
                "grouping": grole, "grouping_column": gcol,
                "label_type": lrole, "label_column": lcol,
                "n_group_levels": counts.shape[0], "n_labels": counts.shape[1],
                "n_cells": int(counts.to_numpy().sum()),
                "cramers_v": cramers_v(counts),
            })

    if long_rows:
        man.write(pd.concat(long_rows, ignore_index=True), "composition_long",
                  "Tidy long-form composition for every grouping x label pair",
                  subdir=SUBDIR)
    if assoc_rows:
        assoc = pd.DataFrame(assoc_rows).sort_values("cramers_v", ascending=False)
        man.write(assoc, "association_summary",
                  "Cramer's V per (grouping, label) pair -- ranks which groupings "
                  "restructure composition most",
                  subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key in cli.selected_datasets(args):
        run(key, args)


if __name__ == "__main__":
    main()
