#!/usr/bin/env python3
"""13 - Do gene lists differ between brain regions within one cell type? (human_dev)

Question: is a list higher in, say, telencephalic than in hindbrain neurons --
consistently in both donor sets -- and how much of that could be age, since
regions were not sampled at the same ages?

Method
  Pseudobulk counts per cell class x region (human_dev; cortex has one
  region) -> log2 TMM-CPM within each class across its regions, genes >= 5
  CPM on average, each gene Z-scored across the class's regions. A set's
  score per region = mean Z of its genes; the statistic for region r is its
  score minus the mean of the class's other regions. Null: random sets drawing
  each member from its expression decile (2,000 sets). v2 and v3 are
  combined per class x region (signed Stouffer, BH, tiered).
  Age check: each region's cell-weighted mean age per stratum (05_confounds
  crosstab), reported as the region's age offset from the class's other
  regions next to every finding; and, per set, the Spearman between region
  effects and region age offsets across all class x region tests.

Inputs (csv_exports/):
  human_dev__<chem>/09_pseudobulk/cell_class_x_region__{pseudobulk_counts,group_summary}.csv
  human_dev__<chem>/05_confounds/crosstab_age_x_region.csv
  human_dev__v2/11_panels/panel_coverage.csv; results/08_coexpression_modules/modules.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "13_regional_differences"
TITLE = "Regional differences in gene lists within cell types (human_dev)"
DATASET = "human_dev"
N_RANDOM = 2000
MIN_GENES = 5
MIN_REGIONS = 3
MIN_CPM = 5.0
MAX_LISTED = 8


def region_ages(n: str) -> pd.Series:
    x = C.csv(n, "05_confounds/crosstab_age_x_region.csv", index_col=0)
    ages = x.index.astype(float).to_numpy()
    w = x.to_numpy(float)
    return pd.Series((w * ages[:, None]).sum(axis=0) / w.sum(axis=0), index=x.columns.astype(str))


def per_stratum(out: C.Output, sets: dict, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_region__pseudobulk_counts.csv",
                 f"{n}/05_confounds/crosstab_age_x_region.csv")
        cnt = C.group_matrix(n, "cell_class_x_region", "pseudobulk_counts")
        gs = C.group_summary(n, "cell_class_x_region").set_index("group")
        rage = region_ages(n)
        by_class: dict[str, list[str]] = {}
        for col in cnt.columns:
            cls, reg = (x.strip() for x in col.rsplit("|", 1))
            by_class.setdefault(cls, []).append(col)
        for cls, cols in sorted(by_class.items()):
            if len(cols) < MIN_REGIONS:
                continue
            lc = C.tmm_log_cpm(cnt[cols])
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            Z = C.zscore_rows(lc.to_numpy(float))
            bins = C.level_bins(lc.mean(axis=1), 10)
            pos = {g: i for i, g in enumerate(lc.index)}
            regions = [c.rsplit("|", 1)[1].strip() for c in cols]
            others = [[j for j in range(len(cols)) if j != i] for i in range(len(cols))]

            def contrast(S: np.ndarray) -> np.ndarray:      # (n_sets x regions) -> same
                return np.stack([S[:, i] - S[:, o].mean(axis=1) for i, o in enumerate(others)], axis=1)

            for name, genes in sets.items():
                idx = np.array([pos[g] for g in genes if g in pos])
                if idx.size < MIN_GENES:
                    continue
                obs = contrast(Z[idx].mean(axis=0)[None, :])[0]
                null = contrast(C.set_mean_rows(Z, C.matched_draws(bins, idx, N_RANDOM, rng)))
                for i, reg in enumerate(regions):
                    eff, p, mu, sd = C.null_effect(obs[i], null[:, i])
                    oth_age = np.mean([rage.get(regions[j], np.nan) for j in others[i]])
                    rows.append({"dataset": DATASET, "chemistry": chem, "cell_class": cls, "region": reg,
                                 "gene_set": name, "n_genes": int(idx.size), "n_regions": len(cols),
                                 "n_cells": int(gs.loc[cols[i], "n_cells"]),
                                 "difference_vs_other_regions": obs[i], "effect_vs_null_sd": eff,
                                 "perm_p": p, "region_mean_age": rage.get(reg, np.nan),
                                 "age_offset_weeks": rage.get(reg, np.nan) - oth_age})
    return pd.DataFrame(rows)


def age_check(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (chem, name), g in per.groupby(["chemistry", "gene_set"]):
        g = g.dropna(subset=["effect_vs_null_sd", "age_offset_weeks"])
        if len(g) >= 8:
            rows.append({"chemistry": chem, "gene_set": name, "n_tests": len(g),
                         "spearman_effect_vs_age_offset": float(
                             g.effect_vs_null_sd.corr(g.age_offset_weeks, method="spearman"))})
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sets = C.analysis_gene_sets(DATASET, groups=("ndd",), modules=True, out=out)
    per = per_stratum(out, sets, rng)
    out.write(per, "region_contrast_per_stratum",
              "Per stratum x class x region x set: set score minus other regions of the class, vs "
              "expression-matched random sets; region age offset")
    comb = C.combine_chemistries(per, ["dataset", "cell_class", "region", "gene_set"],
                                 labels=("higher", "lower"),
                                 carry=("difference_vs_other_regions", "age_offset_weeks", "n_genes"))
    out.write(comb, "region_contrast_combined",
              "Per class x region x set: v2 x v3 combined; tier; region age offset per chemistry")
    ac = age_check(per)
    out.write(ac, "region_effect_vs_age",
              "Per set x stratum: Spearman between region effects and region age offsets (a strong "
              "value means region differences may be age differences)")
    figures(out, comb)

    # ---- findings -----------------------------------------------------------
    f = []
    rep = comb[comb.tier != ""] if not comb.empty else comb
    classes = per.cell_class.nunique()
    f.append(f"**Tested**: {classes} cell classes with >= {MIN_REGIONS} regions, "
             f"{per.gene_set.nunique()} gene sets, {len(comb)} class x region x set tests in both "
             f"chemistries; {len(rep)} tiered ({(rep.tier == 'replicated').sum()} replicated).")
    for kind, label in (("list:", "Gene lists"), ("seed:", "Seed NDD panels"), ("module:", "Modules")):
        r = rep[rep.gene_set.str.startswith(kind)] if len(rep) else rep
        if r.empty:
            if kind != "module:":
                f.append(f"**{label}: no region difference within a class replicates.**")
            continue
        parts = []
        for name, g in r.groupby("gene_set", sort=False):
            g = g.reindex(g.stouffer_z.abs().sort_values(ascending=False).index)
            parts.append(f"{name} ({len(g)} tiered): " + ", ".join(
                f"{x.direction} in {x.cell_class} of {x.region} ({x.effect_v2:+.1f}/{x.effect_v3:+.1f}; "
                f"age {x.age_offset_weeks_v2:+.1f}/{x.age_offset_weeks_v3:+.1f} wk)" for x in g.head(3).itertuples()))
        f.append(f"**{label} that differ between regions within a class** (strongest three per set; "
                 "effect in null SDs v2/v3; region age offset vs the class's other regions): "
                 + " | ".join(parts[:MAX_LISTED * 2])
                 + ("" if len(parts) <= MAX_LISTED * 2 else f" | (+{len(parts) - MAX_LISTED * 2} more sets)") + ".")
    if not ac.empty:
        w = ac.pivot(index="gene_set", columns="chemistry", values="spearman_effect_vs_age_offset").dropna()
        strong = w[(w.abs() >= 0.5).all(axis=1) & (np.sign(w["v2"]) == np.sign(w["v3"]))]
        f.append("**Age check**: "
                 + (f"region effects track region age for {', '.join(strong.index[:MAX_LISTED])} "
                    "(|rho| >= 0.5 in both chemistries) -- read their region differences as possibly age."
                    if len(strong) else "no set's region effects track region age (|rho| >= 0.5 in "
                    "both chemistries), so the differences above are not simply age."))

    out.summary(
        TITLE,
        "Within a cell type, do gene lists and modules differ between brain regions, in both "
        "donor sets, and could the difference be age?",
        [f"Pseudobulk per class x region; log2 TMM-CPM within class; genes >= {MIN_CPM:g} CPM; Z "
         "across the class's regions; set score per region minus mean of the class's other "
         f"regions; {N_RANDOM:,} random sets matched on expression decile; classes with >= "
         f"{MIN_REGIONS} regions.",
         "v2 x v3 signed Stouffer, BH, tiered. Region age offsets from 05_confounds "
         "(cell-weighted mean age of the region in the stratum)."],
        f,
        ["Each class x region pseudobulk pools several donors of varying ages; regions were "
         "sampled at different ages, so a region difference can be partly age (see the offsets "
         "and the age check).",
         "Region labels are dissection labels; 'Forebrain' overlaps Telencephalon and Diencephalon.",
         "Pseudobulks pool cells, not donors, so the per-chemistry p-values treat genes as the "
         "random unit; replication across donor sets is the safeguard."],
        ["B1: class x region x age pseudobulks would compare regions at matched ages directly."])


def figures(out: C.Output, comb: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or comb.empty:
        return
    h = comb[comb.gene_set.str.startswith(("list:", "seed:"))].copy()
    if h.empty:
        return
    h["col"] = h.cell_class + " | " + h.region
    m = h.pivot(index="gene_set", columns="col", values="stouffer_z")
    keep = m.columns[(m.abs() >= 3).any(axis=0)]
    m = m[keep] if len(keep) else m
    fig, ax = plt.subplots(figsize=(0.3 * m.shape[1] + 3, 0.35 * m.shape[0] + 2.5))
    im = ax.imshow(m.to_numpy(float), cmap="RdBu_r", vmin=-6, vmax=6, aspect="auto")
    ax.set_xticks(range(m.shape[1]), m.columns, rotation=75, ha="right", fontsize=6)
    ax.set_yticks(range(m.shape[0]), m.index, fontsize=7)
    ax.set_title("human_dev: gene sets by region within class (combined Z; columns with |Z| >= 3)", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout()
    out.figure(fig, "regional_differences", "Gene sets by region within cell class (combined Z)")
    plt.close(fig)


if __name__ == "__main__":
    main()
