#!/usr/bin/env python3
"""14 - Normalisation factors and diagnostics over the pseudobulk matrices.

CPM assumes every library samples the same underlying population. When a few
genes dominate one group -- routine in brain data, where neurons pour out a
handful of transcripts -- that assumption breaks and every *other* gene looks
depressed. This computes TMM factors, which estimate that composition bias, and
the MA diagnostics that reveal it.

Deliberately NOT computed: RPKM/FPKM. Both datasets are 10x UMI data counted
from the 3' end, so there is no transcript-length bias to divide out; applying
RPKM here would create an artefact rather than remove one (lib/bulk_stats.py
carries the full note).

Outputs (csv_exports/<dataset>/14_normalization/)
  tmm_factors_<grouping>.csv       library size, TMM factor, effective size
  ma_diagnostics_<grouping>.csv    median M vs reference; bias_flag per group
  distribution_<grouping>.csv      per-group expression quantiles (box/violin)
  normalization_summary.csv        one row per grouping: spread of factors
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.bulk_stats import ma_statistics, tmm_normalization_factors
from lib.io_utils import Manifest, log

SCRIPT = "14_normalization_diagnostics"
SUBDIR = "14_normalization"
QUANTILES = [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]


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

    summary = []
    for cpath in sorted(pb.glob("*__pseudobulk_counts.csv")):
        grouping = cpath.name.replace("__pseudobulk_counts.csv", "")
        counts = load_matrix(cpath)
        if counts.shape[1] < 2:
            continue
        log(f"  {grouping}: {counts.shape[0]:,} genes x {counts.shape[1]} groups")

        tmm = tmm_normalization_factors(counts)
        tmm.insert(0, "grouping", grouping)
        man.write(tmm, f"tmm_factors_{grouping}",
                  "TMM scaling factors per group (Robinson & Oshlack): library "
                  "size, factor, and composition-corrected effective size",
                  subdir=SUBDIR)

        ma = ma_statistics(counts)
        if not ma.empty:
            ma.insert(0, "grouping", grouping)
            man.write(ma, f"ma_diagnostics_{grouping}",
                      "MA-plot summary vs the reference library; |median M| > 0.2 "
                      "is flagged as a global shift CPM will not fix",
                      subdir=SUBDIR)

        # Per-group distribution of log2 CPM -- the box/violin input.
        lib = counts.sum(axis=0).replace(0, np.nan)
        logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
        q = logcpm.quantile(QUANTILES).T
        q.columns = [f"log2cpm_q{int(c * 100)}" for c in q.columns]
        q.index.name = "group"
        q = q.reset_index()
        q["mean_log2cpm"] = logcpm.mean().values
        q["sd_log2cpm"] = logcpm.std().values
        q["n_genes_detected"] = (counts > 0).sum().values
        q["library_size"] = counts.sum().values
        q.insert(0, "grouping", grouping)
        man.write(q, f"distribution_{grouping}",
                  "Per-group log2 CPM quantiles, mean/SD and detected-gene count "
                  "-- the box/violin summary for spotting a mis-normalised group",
                  subdir=SUBDIR)

        fac = tmm["tmm_factor"].replace([np.inf, -np.inf], np.nan).dropna()
        n_flag = int(ma["bias_flag"].sum()) if "bias_flag" in ma.columns else 0
        summary.append({
            "grouping": grouping, "n_groups": counts.shape[1],
            "tmm_factor_min": float(fac.min()) if len(fac) else np.nan,
            "tmm_factor_max": float(fac.max()) if len(fac) else np.nan,
            "tmm_factor_spread": float(fac.max() / fac.min()) if len(fac) and fac.min() > 0 else np.nan,
            "n_groups_bias_flagged": n_flag,
            "library_size_min": float(counts.sum().min()),
            "library_size_max": float(counts.sum().max()),
        })
        if n_flag:
            log(f"    {n_flag}/{counts.shape[1]} groups show a global shift (|median M| > 0.2)")

    if summary:
        s = pd.DataFrame(summary)
        man.write(s, "normalization_summary",
                  "Per grouping: spread of TMM factors and how many groups show "
                  "composition bias. A large spread means CPM alone is unsafe.",
                  subdir=SUBDIR)
        worst = s.sort_values("tmm_factor_spread", ascending=False).iloc[0]
        log(f"  widest TMM spread: {worst['grouping']} "
            f"({worst['tmm_factor_spread']:.2f}x)")
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key in cli.selected_datasets(args):
        run(key, args)


if __name__ == "__main__":
    main()
