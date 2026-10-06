#!/usr/bin/env python3
"""21 - Age trends inside one brain region (human_dev): removing the region confound.

Question: human_dev's within-class age trends (03, 06) pool brain regions
whose dissection changes with age, so a "trend" can be a shift in region mix.
Inside one region, which genes and lists still change with age, in both donor
sets -- and does restricting to the telencephalon bring human_dev into line
with cortex (which is telencephalic tissue, from a subset of the same donors)?

Method
  Pseudobulk counts per cell class x region x age (09, from this version on)
  -> log2 TMM-CPM within each class x region across its ages (>= 50 cells per
  age point, >= 5 age points), genes >= 5 CPM. Per gene: Spearman with age,
  exact permutation p; v2 x v3 signed Stouffer (weights sqrt(ages)), BH per
  class x region, tiered. Lists and seed NDD panels: mean rho of their genes
  vs random sets matched on expression decile (as 06), combined and tiered.
  Comparisons: (a) with 03's pooled human_dev trends in the same class --
  correlation of combined Z over genes, and how many pooled replicated trends
  hold within the region; (b) with cortex's 03 trends -- does the region-
  restricted human_dev trend agree better with cortex than the pooled one?

Inputs:
  csv_exports/human_dev__<chem>/09_pseudobulk/cell_class_x_region_x_age__{pseudobulk_counts,group_summary}.csv
  results/03_age_trends_within_cell_class/age_trends_combined.csv (comparisons, optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "21_region_age_trends"
TITLE = "Age trends within one brain region (human_dev): the region confound removed"
DATASET = "human_dev"
GROUPING = "cell_class_x_region_x_age"
MIN_AGES = 5
MIN_CPM = 5.0
N_RANDOM = 2000
MAX_LISTED = 8


def per_stratum(out: C.Output, sets: dict, rng: np.random.Generator):
    genes_rows, set_rows = [], []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        p = C.EXPORTS / n / "09_pseudobulk" / f"{GROUPING}__pseudobulk_counts.csv"
        if not p.exists():
            return None, None
        out.used(f"{n}/09_pseudobulk/{GROUPING}__pseudobulk_counts.csv")
        cnt = C.group_matrix(n, GROUPING, "pseudobulk_counts")
        parts = {}
        for col in cnt.columns:
            cls, reg, age = C.parse_group(col, 3)
            try:
                if C.excluded(DATASET, age=float(age)):
                    continue
                parts.setdefault((cls, reg), []).append((float(age), col))
            except ValueError:
                continue
        for (cls, reg), items in sorted(parts.items()):
            if len(items) < MIN_AGES:
                continue
            items.sort()
            cols = [c for _, c in items]
            ages = np.array([a for a, _ in items])
            lc = C.tmm_log_cpm(cnt[cols])
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            t = C.gene_age_trends(lc, ages)
            t.insert(0, "region", reg)
            t.insert(0, "cell_class", cls)
            t.insert(0, "chemistry", chem)
            t["n_ages"] = len(ages)
            t["ages"] = f"{ages.min():g}-{ages.max():g}"
            genes_rows.append(t)
            st = C.set_shift_test(t.set_index("gene").rho, t.set_index("gene").mean_log2cpm, sets, rng, N_RANDOM)
            if not st.empty:
                st.insert(0, "region", reg)
                st.insert(0, "cell_class", cls)
                st.insert(0, "chemistry", chem)
                st["n_ages"] = len(ages)
                set_rows.append(st)
    if not genes_rows:
        return pd.DataFrame(), pd.DataFrame()
    g = pd.concat(genes_rows, ignore_index=True).assign(dataset=DATASET)
    s = pd.concat(set_rows, ignore_index=True).assign(dataset=DATASET) if set_rows else pd.DataFrame()
    return g, s


def combine_genes(per: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for (cls, reg), g in per.groupby(["cell_class", "region"], sort=False):
        c = C.combine_chemistries(g, ["dataset", "cell_class", "region", "gene"], effect="rho",
                                  weight="n_ages", labels=("rises with age", "falls with age"),
                                  carry=("rho", "ages"))
        parts.append(c)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def comparisons(out: C.Output, comb: pd.DataFrame) -> pd.DataFrame:
    p03 = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_combined.csv"
    if not p03.exists() or comb.empty:
        return pd.DataFrame()
    out.used("results/03_age_trends_within_cell_class/age_trends_combined.csv")
    t03 = pd.read_csv(p03, low_memory=False, dtype={"panels": str})
    pooled = t03[t03.dataset == DATASET].set_index(["cell_class", "gene"])
    cortex = t03[t03.dataset == "cortex"].set_index(["cell_class", "gene"])
    rows = []
    for (cls, reg), g in comb.groupby(["cell_class", "region"]):
        g = g.set_index("gene")
        row = {"cell_class": cls, "region": reg, "n_genes": len(g),
               "n_tiered_within_region": int((g.tier != "").sum())}
        if cls in pooled.index.get_level_values(0):
            pc = pooled.loc[cls]
            shared = g.index.intersection(pc.index)
            row["spearman_z_vs_pooled_human_dev"] = float(g.loc[shared, "stouffer_z"].corr(
                pc.loc[shared, "stouffer_z"], method="spearman")) if len(shared) > 20 else np.nan
            rep = pc[pc.tier == "replicated"].index.intersection(shared)
            same = (np.sign(g.loc[rep, "stouffer_z"]) == np.sign(pc.loc[rep, "stouffer_z"]))
            row["pooled_replicated_genes"] = len(rep)
            row["of_which_same_sign_here"] = int(same.sum())
            row["of_which_tiered_here"] = int((same & (g.loc[rep, "tier"] != "")).sum())
        if cls in cortex.index.get_level_values(0):
            cc = cortex.loc[cls]
            shared = g.index.intersection(cc.index)
            if len(shared) > 20:
                row["spearman_z_vs_cortex"] = float(g.loc[shared, "stouffer_z"].corr(
                    cc.loc[shared, "stouffer_z"], method="spearman"))
                if cls in pooled.index.get_level_values(0):
                    pc = pooled.loc[cls]
                    sh2 = shared.intersection(pc.index)
                    row["spearman_pooled_human_dev_vs_cortex"] = float(pc.loc[sh2, "stouffer_z"].corr(
                        cc.loc[sh2, "stouffer_z"], method="spearman")) if len(sh2) > 20 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sets = C.analysis_gene_sets(DATASET, groups=("ndd",), modules=False, out=out)
    per, sets_per = per_stratum(out, sets, rng)
    if per is None:
        out.summary(TITLE, "Skipped: no class x region x age pseudobulk yet.",
                    [f"Needs 09_pseudobulk/{GROUPING}__pseudobulk_counts.csv."],
                    ["**Not run**: script 09 writes this grouping from this version of the repo on; re-run "
                     "stage 2 (slurm_02_pseudobulk.sh), then this analysis."], [], [])
        return
    if per.empty:
        out.summary(TITLE, "No class x region has enough age points.",
                    [f"Needs >= {MIN_AGES} age points of >= {C.MIN_CELLS} cells per class x region."],
                    [f"**Nothing testable**: no class x region reaches {MIN_AGES} age points in a chemistry."],
                    [], [])
        return
    comb = combine_genes(per)
    scomb = C.combine_chemistries(sets_per, ["dataset", "cell_class", "region", "gene_set"], weight="n_ages",
                                  labels=("rises with age", "falls with age"), carry=("mean",)) \
        if not sets_per.empty else pd.DataFrame()
    cmp = comparisons(out, comb)
    out.write(per, "gene_age_trends_per_stratum", "Per stratum x class x region x gene: Spearman with age")
    out.write(comb, "gene_age_trends_combined", "Per class x region x gene: v2 x v3 combined; BH per class x region; tier")
    out.write(sets_per, "set_age_trends_per_stratum", "Per stratum x class x region x set: mean rho vs matched null")
    out.write(scomb, "set_age_trends_combined", "Per class x region x set: v2 x v3 combined; tier")
    out.write(cmp, "comparison_with_03", "Per class x region: agreement with 03's pooled human_dev and cortex trends")

    f = []
    tested = comb.groupby(["cell_class", "region"]).size() if not comb.empty else pd.Series(dtype=int)
    f.append(f"**Tested**: {len(tested)} class x region combinations with >= {MIN_AGES} age points in both "
             "chemistries: " + ", ".join(f"{c} / {r}" for c, r in tested.index[:20]) + ".")
    if not comb.empty:
        tiers = comb[comb.tier != ""].groupby(["cell_class", "region", "direction"]).size().unstack(fill_value=0)
        f.append("**Genes with a replicated age trend inside one region** (rises / falls): "
                 + ("; ".join(f"{c} / {r}: {int(row.get('rises with age', 0))} / {int(row.get('falls with age', 0))}"
                              for (c, r), row in tiers.iterrows()) if len(tiers) else "none") + ".")
    if not cmp.empty:
        x = cmp.dropna(subset=["spearman_z_vs_cortex", "spearman_pooled_human_dev_vs_cortex"], how="any") \
            if {"spearman_z_vs_cortex", "spearman_pooled_human_dev_vs_cortex"} <= set(cmp.columns) else pd.DataFrame()
        if len(x):
            f.append("**Agreement with cortex** (Spearman of combined Z over genes; region-restricted vs pooled "
                     "human_dev): " + "; ".join(f"{r.cell_class} / {r.region} {r.spearman_z_vs_cortex:+.2f} vs "
                                                f"{r.spearman_pooled_human_dev_vs_cortex:+.2f}" for r in x.itertuples())
                     + ". Higher within the telencephalon means the pooled trends carried region mix.")
        if "pooled_replicated_genes" in cmp:
            y = cmp.dropna(subset=["pooled_replicated_genes"])
            f.append("**03's pooled human_dev replicated trends inside one region** (same sign / tiered of all): "
                     + "; ".join(f"{r.cell_class} / {r.region} {int(r.of_which_same_sign_here)} / "
                                 f"{int(r.of_which_tiered_here)} of {int(r.pooled_replicated_genes)}"
                                 for r in y.itertuples() if r.pooled_replicated_genes > 0) + ".")
    if not scomb.empty:
        rep = scomb[scomb.tier != ""]
        f.append("**Lists and panels moving with age inside one region** (mean rho v2/v3): "
                 + ("; ".join(f"{r.gene_set} {r.direction} in {r.cell_class} / {r.region} ({r.mean_v2:+.2f}/"
                              f"{r.mean_v3:+.2f}; {r.tier})" for r in rep.head(MAX_LISTED * 2).itertuples())
                    if len(rep) else "none replicated") + ".")
    out.summary(
        TITLE,
        "Inside a single brain region, which genes and lists change with age in both donor sets, and "
        "does removing the region confound bring human_dev's trends into line with cortex?",
        [f"Pseudobulk per class x region x age (09); log2 TMM-CPM within class x region; >= {MIN_AGES} ages of "
         f">= {C.MIN_CELLS} cells; genes >= {MIN_CPM:g} CPM; Spearman with age, exact permutation; v2 x v3 "
         "signed Stouffer, BH per class x region, tiered.",
         f"Lists / seed NDD panels: mean rho vs {N_RANDOM:,} random sets matched on expression decile.",
         "Comparison with 03 (pooled human_dev, and cortex) by Spearman of combined Z over shared genes."],
        f,
        ["Each region is sampled at fewer ages than its class overall, so power is lower than in 03.",
         "Region labels are dissection labels; within-region heterogeneity (e.g. dorsal vs ventral "
         "telencephalon) remains."],
        ["Subregion x age pseudobulks for the telencephalon (cortex vs ganglionic eminences)."])


if __name__ == "__main__":
    main()
