#!/usr/bin/env python3
"""23 - Is a gene list active in all cells of a type, or in a subset? (per-cell scores)

Question: a pseudobulk mean cannot tell "every cell expresses the list a
little more" from "a subset expresses it a lot". From per-cell module scores
(stage-2 script 22): how is each list's score spread across the cells of a
class, compared with random lists of matched expression; which share of cells
is clearly active; and does that share change with age or cell-cycle phase?

Method (per stratum; v2 and v3 as replicates)
  Per class, from score histograms: the threshold T = 95th percentile of the
  pooled scores of the list's three matched random programmes in the same
  cells. Active share = share of cells scoring above T (5% expected for a
  list that behaves like random genes). Spread ratio = SD of the list's score
  / mean SD of its random programmes (> 1: cells differ more than random
  genes would make them differ -- a subset, or a gradient). Shift = mean
  score minus the random programmes' mean.
  Pattern: 'subset' = spread ratio >= 1.5 in both chemistries; 'broad' = a
  positive shift and an active share of at least 10% (twice what random
  lists reach) with spread ratio < 1.5, in both; otherwise unresolved.
  Age: active share per class x age (threshold per group), Spearman with age,
  exact permutation, v2 x v3 combined, tiered. Phase: active share per class
  x phase.

Inputs (csv_exports/<ds>__<chem>/22_cell_programs/):
  programs.csv, program_score_summary.csv, program_score_hist.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "23_cell_level_activity"
TITLE = "Gene lists at single-cell level: broad or subset activity, by age and phase"
TAIL = 0.95
SUBSET_SPREAD = 1.5
BROAD_SHARE = 2 * (1 - TAIL)   # 'broad' needs at least twice the share random lists reach
MIN_GROUP_CELLS = 50
MIN_AGES = 5
MAX_LISTED = 8


def load(n: str, out: C.Output):
    base = C.EXPORTS / n / "22_cell_programs"
    if not (base / "program_score_hist.csv").exists():
        return None
    for f in ("programs", "program_score_summary", "program_score_hist"):
        out.used(f"{n}/22_cell_programs/{f}.csv")
    progs = pd.read_csv(base / "programs.csv")
    summ = pd.read_csv(base / "program_score_summary.csv")
    hist = pd.read_csv(base / "program_score_hist.csv")
    H = hist.pivot_table(index=["grouping", "group", "program"], columns="bin_lo", values="count",
                         aggfunc="sum", fill_value=0)
    H = H.reindex(sorted(H.columns), axis=1)
    return progs, summ, H


def activity(progs: pd.DataFrame, summ: pd.DataFrame, H: pd.DataFrame, grouping: str) -> pd.DataFrame:
    lo = np.array(H.columns, dtype=float)
    width = float(np.nanmin(np.diff(lo[np.isfinite(lo)])))
    hi = np.where(np.isfinite(lo), lo + width, np.nanmin(lo[np.isfinite(lo)]))
    s = summ[summ.grouping == grouping].set_index(["group", "program"])
    real = progs[progs.kind == "programme"].program
    rand = progs[progs.kind == "random"].groupby("source").program.apply(list)
    rows = []
    groups = H.loc[grouping].index.get_level_values("group").unique() if grouping in H.index.get_level_values(0) else []
    for grp in groups:
        Hg = H.loc[(grouping, grp)]
        for p in real:
            if p not in Hg.index or p not in rand.index:
                continue
            rr = [r for r in rand[p] if r in Hg.index]
            if not rr:
                continue
            hp = Hg.loc[p].to_numpy(float)
            hr = Hg.loc[rr].to_numpy(float).sum(axis=0)
            n = hp.sum()
            if n < MIN_GROUP_CELLS or hr.sum() == 0:
                continue
            k = int(np.searchsorted(np.cumsum(hr) / hr.sum(), TAIL))
            thr = hi[min(k, len(hi) - 1)]
            active = float(hp[lo >= thr].sum() / n)
            sd_p = float(s.loc[(grp, p), "sd"])
            sd_r = float(np.mean([s.loc[(grp, r), "sd"] for r in rr]))
            rows.append({"group": grp, "program": p, "n_cells": int(n), "threshold": thr,
                         "active_share": active, "spread_ratio": sd_p / sd_r if sd_r > 0 else np.nan,
                         "shift": float(s.loc[(grp, p), "mean"] - np.mean([s.loc[(grp, r), "mean"] for r in rr]))})
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    cls_rows, age_rows, phase_rows = [], [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        loaded = load(n, out)
        if loaded is None:
            continue
        progs, summ, H = loaded
        a = activity(progs, summ, H, "cell_class")
        if len(a):
            cls_rows.append(a.rename(columns={"group": "cell_class"}).assign(dataset=ds, chemistry=chem))
        b = activity(progs, summ, H, "cell_class_x_age")
        if len(b):
            parts = b.group.map(lambda g: C.parse_group(g, 2))
            b["cell_class"], b["age"] = parts.str[0], pd.to_numeric(parts.str[1], errors="coerce")
            b = b[~b.age.map(lambda x: C.excluded(ds, age=x) if pd.notna(x) else True)]
            for (cls, p), g in b.groupby(["cell_class", "program"]):
                if g.age.nunique() < MIN_AGES:
                    continue
                g = g.sort_values("age")
                rho = float(C.spearman_rows(g.active_share.to_numpy()[None, :], g.age.to_numpy())[0])
                pp, _ = C.spearman_perm_p(np.array([rho]), g.age.to_numpy())
                age_rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "program": p,
                                 "n_ages": len(g), "rho_active_share_vs_age": rho, "perm_p": float(pp[0]),
                                 "active_share_youngest": float(g.active_share.iloc[0]),
                                 "active_share_oldest": float(g.active_share.iloc[-1])})
        c = activity(progs, summ, H, "cell_class_x_phase")
        if len(c):
            parts = c.group.map(lambda g: C.parse_group(g, 2))
            c["cell_class"], c["phase"] = parts.str[0], parts.str[1]
            phase_rows.append(c.assign(dataset=ds, chemistry=chem))
    if not cls_rows:
        out.summary(TITLE, "Skipped: no per-cell programme scores yet.",
                    ["Needs 22_cell_programs/ from stage 2 (script 22)."],
                    ["**Not run**: script 22 is new; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis."],
                    [], [])
        return
    cls = pd.concat(cls_rows, ignore_index=True)
    w = cls.pivot_table(index=["dataset", "cell_class", "program"], columns="chemistry",
                        values=["active_share", "spread_ratio", "shift"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    if {"spread_ratio_v2", "spread_ratio_v3"} <= set(w.columns):
        both = lambda col, cond: cond(w[f"{col}_v2"]) & cond(w[f"{col}_v3"])
        w["pattern"] = np.where(both("spread_ratio", lambda x: x >= SUBSET_SPREAD), "subset",
                                np.where(both("active_share", lambda x: x >= BROAD_SHARE)
                                         & both("shift", lambda x: x > 0)
                                         & both("spread_ratio", lambda x: x < SUBSET_SPREAD), "broad", ""))
    ages = pd.DataFrame(age_rows)
    acomb = C.combine_chemistries(ages, ["dataset", "cell_class", "program"], effect="rho_active_share_vs_age",
                                  weight="n_ages", labels=("more cells active with age", "fewer cells active with age"),
                                  carry=("active_share_youngest", "active_share_oldest")) if len(ages) else pd.DataFrame()
    phase = pd.concat(phase_rows, ignore_index=True) if phase_rows else pd.DataFrame()
    out.write(cls, "activity_per_stratum", "Per stratum x class x list: active share (above the random lists' 95th "
              "percentile), spread ratio, shift")
    out.write(w, "activity_by_class", "Per class x list: metrics per chemistry and pattern (subset / broad)")
    out.write(ages, "active_share_vs_age_per_stratum", "Per stratum x class x list: active share vs age")
    out.write(acomb, "active_share_vs_age_combined", "v2 x v3 combined; tier")
    out.write(phase, "activity_by_phase", "Per stratum x class x phase x list: active share, spread ratio")

    f = []
    for ds in C.DATASETS:
        g = w[(w.dataset == ds) & w.program.str.startswith("list:")] if "pattern" in w else w.iloc[0:0]
        if g.empty:
            g = w[(w.dataset == ds)] if "pattern" in w else g
        sub = g[g.pattern == "subset"] if "pattern" in g else g.iloc[0:0]
        broad = g[g.pattern == "broad"] if "pattern" in g else g.iloc[0:0]
        f.append(f"**{ds}: lists active in a subset of a class's cells** (spread ratio v2/v3; active share v2/v3): "
                 + ("; ".join(f"{r.program} in {r.cell_class} ({r.spread_ratio_v2:.1f}/{r.spread_ratio_v3:.1f}; "
                              f"{r.active_share_v2:.0%}/{r.active_share_v3:.0%})"
                              for r in sub.sort_values("spread_ratio_v2", ascending=False).head(MAX_LISTED).itertuples())
                    if len(sub) else "none") + ".")
        f.append(f"**{ds}: lists raised broadly across a class** (active share v2/v3): "
                 + ("; ".join(f"{r.program} in {r.cell_class} ({r.active_share_v2:.0%}/{r.active_share_v3:.0%})"
                              for r in broad.sort_values("active_share_v2", ascending=False).head(MAX_LISTED).itertuples())
                    if len(broad) else "none") + ".")
    ctrl = w[w.program.isin(["seed:g2m_phase", "seed:s_phase"])] if "pattern" in w else w.iloc[0:0]
    if len(ctrl):
        f.append("**Positive control** (phase panels should be 'subset' in progenitors -- only cycling cells "
                 "express them): " + "; ".join(f"{r.dataset} {r.program} in {r.cell_class}: {r.pattern or 'unresolved'} "
                                               f"(spread {r.spread_ratio_v2:.1f}/{r.spread_ratio_v3:.1f})"
                                               for r in ctrl[ctrl.cell_class.isin(["Radial glia", "Neuronal IPC"])].itertuples()) + ".")
    if not acomb.empty:
        rep = acomb[(acomb.tier != "") & acomb.program.str.startswith(("list:", "seed:"))]
        f.append("**Share of active cells changing with age** (youngest -> oldest, v2 / v3): "
                 + ("; ".join(f"{r.program} {r.direction} in {r.dataset} {r.cell_class} "
                              f"({r.active_share_youngest_v2:.0%}->{r.active_share_oldest_v2:.0%} / "
                              f"{r.active_share_youngest_v3:.0%}->{r.active_share_oldest_v3:.0%}; {r.tier})"
                              for r in rep.head(MAX_LISTED * 2).itertuples()) if len(rep) else "none replicated") + ".")
    out.summary(
        TITLE,
        "For each gene list: is it active in all cells of a class or in a subset, compared with random "
        "lists of matched expression; and does the share of active cells change with age or cell-cycle phase?",
        ["Per-cell module scores from stage 2 (script 22): mean log1p(CP10K) of the list's genes minus "
         "expression-matched control genes; three matched random programmes per list as reference.",
         f"Active share = cells above the random programmes' {TAIL:.0%} quantile in the same group; spread "
         f"ratio = list SD / random SD; 'subset' = spread ratio >= {SUBSET_SPREAD} in both chemistries.",
         f"Age: active share per class x age (>= {MIN_GROUP_CELLS} cells), Spearman with age, exact "
         "permutation, v2 x v3 signed Stouffer, tiered."],
        f,
        ["Module scores average over a list; a subset pattern can mean a sub-type expressing the list or a "
         "gradient across the class.",
         "Sparse single-cell counts make every score noisy; the random programmes share that noise, which is "
         "why they are the reference. They are matched on total UMIs over the whole dataset, not within the "
         "class, so a list of genes typical for the class can sit somewhat above them; hence 'broad' needs "
         "twice the random share."],
        ["Link subset-active cells to clusters or pseudotime (25) to name the subset."])


if __name__ == "__main__":
    main()
