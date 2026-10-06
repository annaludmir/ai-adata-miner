#!/usr/bin/env python3
"""28 - Do the main results hold with stricter cell QC?

Question: the files are already QC-filtered, but tails remain (high
mitochondrial or unspliced fraction, few genes, likely doublets). Rerunning
the pipeline with stricter thresholds (exclusions_strict.csv; AIM_STRICT=true)
gives a second set of results. Which tiered results survive, which are lost,
and do effect sizes agree?

Method
  For each key table, the standard results (this run) are joined with the
  strict-QC results (AIM_STRICT_RESULTS) on the table's keys. Reported: tiered
  results in each run, kept in both, lost and gained under strict QC; among
  results tiered in the standard run, the share with the same sign under strict
  QC; and the Spearman correlation of the combined statistic (Stouffer Z, or
  the summed per-chemistry effects) over all shared rows.

Inputs: results/<analysis>/<table>.csv in both result folders. Skipped inside a
strict run itself, or when no strict results exist.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "28_strict_qc_comparison"
TITLE = "Robustness to stricter cell QC: standard vs strict-QC results"
STRICT = Path(os.environ.get("AIM_STRICT_RESULTS", "/miridan-data/annaludmir/aim_downstream_strict/results"))
SPECS = [
    ("02 composition vs age", "02_composition_vs_age/composition_trends_replicated.csv", ["dataset", "scope", "cell_class"]),
    ("03 gene age trends", "03_age_trends_within_cell_class/age_trends_combined.csv", ["dataset", "cell_class", "gene"]),
    ("06 list cell-class preference", "06_gene_list_landscape/list_class_preference_combined.csv",
     ["dataset", "gene_list", "cell_class"]),
    ("06 list age coordination", "06_gene_list_landscape/list_age_coordination_combined.csv",
     ["dataset", "cell_class", "gene_list"]),
    ("07 list coherence", "07_gene_list_coherence/coherence_combined.csv", ["dataset", "context", "gene_list"]),
    ("09 proliferation trajectories", "09_cell_cycle_programs/proliferation_trajectories_combined.csv",
     ["dataset", "cell_class", "metric"]),
    ("10 splicing induction", "10_splicing_dynamics/set_induction_by_cell_class.csv", ["cell_class", "gene_set"]),
    ("12 sub-type abundance", "12_cluster_abundance_vs_age/cluster_abundance_combined.csv", ["dataset", "parent", "cluster"]),
]
MAX_LISTED = 6


def stat(df: pd.DataFrame) -> pd.Series:
    if "stouffer_z" in df:
        return df["stouffer_z"]
    if {"effect_v2", "effect_v3"} <= set(df.columns):
        return df["effect_v2"] + df["effect_v3"]
    return pd.Series(np.nan, index=df.index)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    skip = None
    if os.environ.get("AIM_STRICT", "false") == "true":
        skip = "this is the strict-QC run itself"
    elif not STRICT.exists() or STRICT.resolve() == C.RESULTS.resolve():
        skip = f"no strict-QC results at {STRICT}"
    if skip:
        out.summary(TITLE, f"Skipped: {skip}.",
                    ["Compares this run's results with a strict-QC rerun (AIM_STRICT_RESULTS)."],
                    [f"**Not run**: {skip}. To produce them: `AIM_STRICT=true ./submit_all.sh` (exclusions_strict.csv; "
                     "exports to csv_exports_strict/, results to aim_downstream_strict/), then rerun stage 4 here."],
                    [], [])
        return
    rows, lost_rows = [], []
    for label, rel, keys in SPECS:
        a, b = C.RESULTS / rel, STRICT / rel
        if not (a.exists() and b.exists()):
            continue
        out.used(f"results/{rel}", f"strict/{rel}")
        s = pd.read_csv(a, low_memory=False)
        t = pd.read_csv(b, low_memory=False)
        if not set(keys) <= set(s.columns) or not set(keys) <= set(t.columns):
            continue
        s, t = s.assign(_stat=stat(s), tier=s.get("tier", pd.Series("", index=s.index)).fillna("")), \
            t.assign(_stat=stat(t), tier=t.get("tier", pd.Series("", index=t.index)).fillna(""))
        m = s[keys + ["_stat", "tier"]].merge(t[keys + ["_stat", "tier"]], on=keys, how="outer",
                                               suffixes=("_standard", "_strict"))
        m[["tier_standard", "tier_strict"]] = m[["tier_standard", "tier_strict"]].fillna("")
        ts, tt = m.tier_standard != "", m.tier_strict != ""
        both = m.dropna(subset=["_stat_standard", "_stat_strict"])
        same = (np.sign(m.loc[ts, "_stat_standard"]) == np.sign(m.loc[ts, "_stat_strict"]))
        rows.append({"table": label, "n_rows_shared": len(both), "tiered_standard": int(ts.sum()),
                     "tiered_strict": int(tt.sum()), "tiered_both": int((ts & tt).sum()),
                     "lost_under_strict": int((ts & ~tt).sum()), "gained_under_strict": int((~ts & tt).sum()),
                     "same_sign_among_standard_tiered": float(same.mean()) if len(same) else np.nan,
                     "spearman_stat": float(both._stat_standard.corr(both._stat_strict, method="spearman"))
                     if len(both) > 10 else np.nan})
        lost = m[ts & ~tt].copy()
        if len(lost):
            lost.insert(0, "table", label)
            lost["key"] = lost[keys].astype(str).agg(" | ".join, axis=1)
            lost_rows.append(lost[["table", "key", "_stat_standard", "_stat_strict", "tier_standard"]])
    res = pd.DataFrame(rows)
    lost = pd.concat(lost_rows, ignore_index=True) if lost_rows else pd.DataFrame()
    out.write(res, "strict_qc_agreement", "Per table: tiered results in each run, kept / lost / gained, sign "
              "agreement and correlation of the combined statistic")
    out.write(lost, "lost_under_strict_qc", "Results tiered with standard QC but not with strict QC")
    f = []
    if res.empty:
        f.append("**No comparable tables found** in both result folders.")
    for r in res.itertuples():
        f.append(f"**{r.table}**: {r.tiered_standard} tiered (standard) vs {r.tiered_strict} (strict); kept "
                 f"{r.tiered_both}, lost {r.lost_under_strict}, gained {r.gained_under_strict}; same sign "
                 f"{r.same_sign_among_standard_tiered:.0%} of standard-tiered; Spearman of statistic "
                 f"{r.spearman_stat:+.2f}.")
    if len(lost):
        f.append("**Examples lost under strict QC**: " + "; ".join(
            f"{t}: " + ", ".join(g.key.head(MAX_LISTED)) for t, g in lost.groupby("table", sort=False)) + ".")
    out.summary(
        TITLE,
        "Which tiered results survive stricter cell QC, and do effect sizes agree?",
        [f"Strict QC: exclusions_strict.csv (mitochondrial fraction > 0.02, < 1,100 genes, unspliced fraction "
         f"> 0.37, cortex doublet score > 0.17 -- roughly each metric's 5% tail). Strict results from {STRICT}.",
         "Per table: join on keys; tiers and combined statistic (Stouffer Z or summed effects) compared."],
        f,
        ["Strict QC removes cells, so it also lowers power: a lost result can be a weaker test, not a QC "
         "artefact; check whether its statistic kept its sign.",
         "Thresholds are set from the exported QC quantiles; they are a robustness probe, not a recommendation."],
        ["Re-derive thresholds per cell class (neurons and progenitors differ in genes per cell)."])


if __name__ == "__main__":
    main()
