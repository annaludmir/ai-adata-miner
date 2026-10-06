#!/usr/bin/env python3
"""24 - Do a list's genes co-vary from cell to cell inside one cell type?

Question: 07 found which lists co-express across clusters -- partly because
clusters differ in identity. Inside one cell class, across single cells, do
the list's genes still rise and fall together (a programme switched on and off
in individual cells), consistently in both donor sets? Which members carry it?

Method
  From stage-2 script 22: per class, correlations of log1p(CP10K) across up to
  20,000 cells with sequencing depth partialled out; a list's coherence is the
  mean pairwise correlation of its genes, against 500 random gene sets matched
  on mean expression in that class. Calibration first: script 22 tests random
  programmes the same way, and their p-values say whether the null holds on
  these data. Each list's effect is also ranked against the random programmes'
  effects in the same class (empirical p), which stays valid even if the null
  is off. v2 x v3: signed Stouffer of the matched-null test, BH per dataset,
  tiered -- and a tier also needs empirical p <= 0.05 in both chemistries.
  Hubs: members with the highest within-cell connectivity in both chemistries.
  Comparison with 07's across-cluster coherence.

Inputs:
  csv_exports/<ds>__<chem>/22_cell_programs/within_cell_{coherence,connectivity}.csv
  results/07_gene_list_coherence/coherence_combined.csv (comparison, optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "24_within_cell_coexpression"
TITLE = "Within-cell co-expression of gene lists inside cell types"
N_HUBS = 5
MAX_LISTED = 8


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    coh, con = [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        base = C.EXPORTS / n / "22_cell_programs"
        if not (base / "within_cell_coherence.csv").exists():
            continue
        out.used(f"{n}/22_cell_programs/within_cell_coherence.csv", f"{n}/22_cell_programs/within_cell_connectivity.csv")
        c = pd.read_csv(base / "within_cell_coherence.csv").assign(dataset=ds, chemistry=chem)
        if "kind" not in c:
            c["kind"] = "programme"
        coh.append(c)
        if (base / "within_cell_connectivity.csv").exists():
            con.append(pd.read_csv(base / "within_cell_connectivity.csv").assign(dataset=ds, chemistry=chem))
    if not coh:
        out.summary(TITLE, "Skipped: no within-cell co-expression yet.",
                    ["Needs 22_cell_programs/within_cell_coherence.csv from stage 2 (script 22)."],
                    ["**Not run**: script 22 is new; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis."],
                    [], [])
        return
    coh = pd.concat(coh, ignore_index=True)
    con = pd.concat(con, ignore_index=True) if con else pd.DataFrame()

    rand = coh[coh.kind == "random"]
    calib = (rand.groupby(["dataset", "chemistry"]).perm_p
             .agg(n="size", frac_p05=lambda p: (p < 0.05).mean(), frac_p01=lambda p: (p < 0.01).mean())
             .reset_index()) if len(rand) else pd.DataFrame()
    real = coh[coh.kind == "programme"].copy()
    emp = []
    for r in real.itertuples():
        ref = rand[(rand.dataset == r.dataset) & (rand.chemistry == r.chemistry) & (rand.cell_class == r.cell_class)]
        e = ref.effect_vs_null_sd.dropna()
        emp.append((np.sum(e >= r.effect_vs_null_sd) + 1) / (len(e) + 1) if len(e) else np.nan)
    real["empirical_p_vs_random_programmes"] = emp
    comb = C.combine_chemistries(real, ["dataset", "cell_class", "program"], labels=("coherent", "anti-coherent"),
                                 carry=("coherence", "n_genes", "empirical_p_vs_random_programmes"))
    if not comb.empty:
        # a tier also needs the list to beat the random programmes in both chemistries,
        # which stays valid when the matched null is miscalibrated
        beats = (comb.empirical_p_vs_random_programmes_v2 <= 0.05) & (comb.empirical_p_vs_random_programmes_v3 <= 0.05)
        comb["tier"] = np.where(beats, comb.tier, "")
    hubs = pd.DataFrame()
    if len(con):
        w = con.pivot_table(index=["dataset", "cell_class", "program", "gene"], columns="chemistry",
                            values="connectivity").dropna().reset_index()
        if {"v2", "v3"} <= set(w.columns):
            w["mean_connectivity"] = (w.v2 + w.v3) / 2
            hubs = (w.sort_values("mean_connectivity", ascending=False)
                    .groupby(["dataset", "cell_class", "program"]).head(N_HUBS))
    cmp07 = pd.DataFrame()
    p07 = C.RESULTS / "07_gene_list_coherence" / "coherence_combined.csv"
    if p07.exists() and not comb.empty:
        out.used("results/07_gene_list_coherence/coherence_combined.csv")
        c07 = pd.read_csv(p07)
        c07 = c07[c07.context == "within_class"].assign(program=lambda d: "list:" + d.gene_list)
        best = comb.sort_values("stouffer_z", ascending=False).drop_duplicates(["dataset", "program"])
        cmp07 = best.merge(c07[["dataset", "program", "tier", "effect_v2", "effect_v3"]],
                           on=["dataset", "program"], suffixes=("_within_cell", "_07_within_class"))
    out.write(coh, "within_cell_coherence_per_stratum", "Per stratum x class x programme (lists, panels and the "
              "random programmes used for calibration)")
    out.write(calib, "null_calibration", "Share of random programmes with p < 0.05 / 0.01 (should be ~0.05 / 0.01)")
    out.write(comb, "within_cell_coherence_combined", "Per class x list: v2 x v3 combined; tier; empirical p vs random programmes")
    out.write(hubs, "within_cell_hubs", f"Per class x list: the {N_HUBS} members most connected to the rest, mean over chemistries")
    out.write(cmp07, "comparison_with_07", "Strongest within-cell result per list vs 07's within-class (across-cluster) coherence")

    f = []
    if len(calib):
        f.append("**Null calibration** (random programmes with p < 0.05 / < 0.01; nominal 5% / 1%): "
                 + "; ".join(f"{r.dataset} {r.chemistry} {r.frac_p05:.1%} / {r.frac_p01:.1%} (n = {r.n})"
                             for r in calib.itertuples())
                 + (". Above nominal, so lean on the empirical p against random programmes."
                    if (calib.frac_p05 > 0.08).any() else ". Close to nominal."))
    for ds in C.DATASETS:
        g = comb[(comb.dataset == ds) & (comb.tier != "") & (comb.direction == "coherent")] if not comb.empty else comb
        f.append(f"**{ds}: lists whose genes co-vary within single cells** (coherence v2/v3; empirical p vs random "
                 "programmes v2/v3): " + ("; ".join(
                     f"{r.program} in {r.cell_class} ({r.coherence_v2:.3f}/{r.coherence_v3:.3f}; "
                     f"{C.fmt_p(r.empirical_p_vs_random_programmes_v2)}/{C.fmt_p(r.empirical_p_vs_random_programmes_v3)}"
                     f"; {r.tier})" for r in g.head(MAX_LISTED * 2).itertuples()) if len(g) else "none replicated") + ".")
    if len(hubs):
        lists = hubs[hubs.program.str.startswith("list:")]
        tiered = set(map(tuple, comb.loc[comb.tier != "", ["dataset", "cell_class", "program"]].to_numpy()))
        sel = lists[[tuple(x) in tiered for x in lists[["dataset", "cell_class", "program"]].to_numpy()]]
        if len(sel):
            f.append("**Hub genes of co-varying lists** (highest mean within-cell connectivity): " + "; ".join(
                f"{p} in {ds} {cls}: " + ", ".join(g.gene) for (ds, cls, p), g in
                list(sel.groupby(["dataset", "cell_class", "program"], sort=False))[:MAX_LISTED]) + ".")
    if len(cmp07):
        f.append("**Within cells vs across clusters (07, within class)**: " + "; ".join(
            f"{r.dataset} {r.program}: within-cell {r.tier_within_cell or 'n.s.'} (best in {r.cell_class}), "
            f"07 {r.tier_07_within_class if isinstance(r.tier_07_within_class, str) and r.tier_07_within_class else 'n.s.'}"
            for r in cmp07.head(MAX_LISTED * 2).itertuples()) + ".")
    out.summary(
        TITLE,
        "Inside one cell class, do a list's genes rise and fall together from cell to cell, in both donor "
        "sets; which members carry it; and how does that compare with co-expression across clusters (07)?",
        ["Stage-2 script 22: per class, correlations across up to 20,000 single cells of log1p(CP10K), "
         "sequencing depth partialled out; coherence = mean pairwise correlation; 500 random sets matched on "
         "class mean as null.",
         "Calibration from random programmes tested the same way; each list's effect also ranked against them.",
         "v2 x v3 signed Stouffer, BH per dataset, tiered; hubs by mean connectivity over chemistries."],
        f,
        ["Single-cell correlations are low in absolute terms (sparse counts); compare with the null, not with 07.",
         "Within a class, sub-types still differ; co-variation across cells includes sub-type differences "
         "inside the class, not only switching within a cell."],
        ["Repeat within single clusters to remove remaining sub-type structure."])


if __name__ == "__main__":
    main()
