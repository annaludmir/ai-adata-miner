#!/usr/bin/env python3
"""15 - Are NDD genes held tighter between donors than comparable genes?

Question: many NDD genes are dosage-sensitive -- losing one copy causes
disease. If their expression is under tight control, it should vary less
between people than that of comparable genes at the same age. Within a cell
class, after removing the age trend, do NDD lists vary less (or more) from
donor to donor than expression-matched genes, in both donor sets?

Method
  Per stratum x cell class with >= 5 age points (an age point is ~one donor):
  log2 TMM-CPM per class x age (as 03), genes >= 5 CPM. Each gene's residual
  SD after a linear fit on age measures its between-donor variability.
  Variability rises steeply as expression falls (counting noise), so each
  gene's log residual SD is turned into a percentile within its expression
  bin (20 quantile bins): 0 = least variable among genes of its level.
  A set's statistic is its genes' mean percentile, against random sets
  drawing each member from the same bin (null mean 0.5). v2 x v3 combined
  (signed Stouffer, weights sqrt(age points)), BH per dataset, tiered.
  Gene table: each gene's percentile per class and stratum, and the genes
  consistently in the tightest fifth in both chemistries.

Inputs:
  csv_exports/<ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,detection_fraction,group_summary}.csv
  csv_exports/<ds>__v2/11_panels/panel_coverage.csv; gene lists; results/08 modules
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "15_donor_variability"
TITLE = "Between-donor variability of NDD genes, against expression-matched genes"
N_RANDOM = 5000
MIN_GENES = 5
MIN_AGES = 5
MIN_CPM = 5.0
N_LEVEL_BINS = 20
TIGHT = 0.2
MAX_LISTED = 8


def variability(lc: pd.DataFrame, ages: np.ndarray) -> pd.DataFrame:
    """Per gene: mean log2 CPM, residual SD after a linear age fit, percentile within level bin."""
    Y = lc.to_numpy(float)
    D = np.stack([np.ones_like(ages, dtype=float), ages.astype(float)], axis=1)
    beta = np.linalg.lstsq(D, Y.T, rcond=None)[0]
    E = Y - C.dot(D, beta).T
    sd = np.sqrt((E ** 2).sum(axis=1) / (len(ages) - 2))
    level = Y.mean(axis=1)
    bins = C.level_bins(level, N_LEVEL_BINS)
    logsd = np.log(np.maximum(sd, 1e-6))
    pct = pd.Series(logsd).groupby(bins).rank(pct=True).to_numpy()
    return pd.DataFrame({"gene": lc.index, "mean_log2cpm": level, "residual_sd": sd,
                         "level_bin": bins, "variability_percentile": pct})


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sets = {ds: C.analysis_gene_sets(ds, groups=("ndd",), modules=True, out=out) for ds in C.DATASETS}
    set_rows, gene_rows = [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv")
        for cls, (lc, _, ages, _) in C.class_age_logcpm(n).items():
            if len(ages) < MIN_AGES:
                continue
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            v = variability(lc, ages)
            gene_rows.append(v.assign(dataset=ds, chemistry=chem, cell_class=cls, n_ages=len(ages)))
            pct, bins = v.variability_percentile.to_numpy(), v.level_bin.to_numpy()
            pos = {g: i for i, g in enumerate(v.gene)}
            for name, genes in sets[ds].items():
                idx = np.array([pos[g] for g in genes if g in pos])
                if idx.size < MIN_GENES:
                    continue
                obs = float(pct[idx].mean())
                null = C.set_mean_rows(pct[:, None], C.matched_draws(bins, idx, N_RANDOM, rng))[:, 0]
                eff, p, mu, sd = C.null_effect(obs, null)
                set_rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "gene_set": name,
                                 "n_genes": int(idx.size), "n_ages": len(ages),
                                 "mean_percentile": obs, "null_mean": mu, "effect_vs_null_sd": eff,
                                 "perm_p": p, "frac_in_tightest_fifth": float((pct[idx] < TIGHT).mean())})
    per = pd.DataFrame(set_rows)
    genes = pd.concat(gene_rows, ignore_index=True)
    comb = C.combine_chemistries(per, ["dataset", "cell_class", "gene_set"], weight="n_ages",
                                 labels=("more variable", "less variable"),
                                 carry=("mean_percentile", "n_genes"))
    out.write(per, "set_variability_per_stratum",
              "Per stratum x class x set: mean variability percentile of its genes vs matched random sets")
    out.write(comb, "set_variability_combined", "v2 x v3 combined; tier; direction vs matched genes")
    w = genes.pivot_table(index=["dataset", "cell_class", "gene"], columns="chemistry",
                          values="variability_percentile").dropna()
    tight = w[(w < TIGHT).all(axis=1)].reset_index()
    out.write(genes, "gene_variability", "Per gene x class x stratum: residual SD after age, percentile "
              "within expression bin (0 = least variable)")
    out.write(tight, "consistently_tight_genes",
              f"Genes in the tightest {TIGHT:.0%} of their expression bin in both chemistries, per class")
    figures(out, comb)

    # ---- findings -----------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        g = comb[comb.dataset == ds] if not comb.empty else comb
        rep = g[g.tier != ""] if len(g) else g
        if rep.empty:
            f.append(f"**{ds}: no set's between-donor variability differs from matched genes in both "
                     "donor sets.**")
            continue
        parts = []
        for direction in ("less variable", "more variable"):
            r = rep[rep.direction == direction]
            if len(r):
                parts.append(f"{direction}: " + "; ".join(
                    f"{x.gene_set} in {x.cell_class} (mean percentile {x.mean_percentile_v2:.2f}/"
                    f"{x.mean_percentile_v3:.2f}; {x.tier})" for x in r.head(MAX_LISTED).itertuples())
                    + ("" if len(r) <= MAX_LISTED else f" (+{len(r) - MAX_LISTED} more)"))
        f.append(f"**{ds}: sets whose between-donor variability differs from matched genes** "
                 "(0.5 = like matched genes): " + " | ".join(parts) + ".")
    f.append("**Consistently tight genes** (tightest fifth of their expression bin in both chemistries): "
             + "; ".join(f"{ds} {cls} {len(g)}" for (ds, cls), g in tight.groupby(["dataset", "cell_class"]))
             + ". See consistently_tight_genes.csv.")
    out.summary(
        TITLE,
        "After removing the age trend, do NDD gene lists vary less (or more) between donors than "
        "genes of the same expression level, within a cell class and in both donor sets?",
        [f"Per stratum x class with >= {MIN_AGES} age points: log2 TMM-CPM per class x age; genes >= "
         f"{MIN_CPM:g} CPM; residual SD after a linear age fit; log SD as a percentile within "
         f"{N_LEVEL_BINS} expression bins.",
         f"Set statistic = mean percentile vs {N_RANDOM:,} random sets matched bin for bin; v2 x v3 "
         "signed Stouffer (weights sqrt(age points)), BH per dataset, tiered."],
        f,
        ["Residual SD mixes biological between-donor variation with dissection and sampling "
         "differences, and with any non-linear age trend.",
         "Co-regulated programmes (cell cycle, neuronal genes, most modules) move together when a "
         "donor's pseudobulk holds more cycling or more mature cells, so they read as 'more "
         "variable'; that is sub-type mix varying between donors, not loose control of each gene. "
         "Genes are treated as independent in the null, which makes such sets reach significance "
         "easily.",
         "An age point is about one donor, so residuals are between donors only approximately "
         "(two donors of one age are pooled).",
         "Dosage sensitivity is one reason for low variability; housekeeping-like stable expression "
         "is another."],
        ["C3: compare with gnomAD LOEUF -- do the tightest genes concentrate among loss-of-function "
         "intolerant genes?"])


def figures(out: C.Output, comb: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or comb.empty:
        return
    for ds in C.DATASETS:
        h = comb[comb.dataset == ds]
        if h.empty:
            continue
        m = h.assign(mp=(h.mean_percentile_v2 + h.mean_percentile_v3) / 2).pivot(
            index="gene_set", columns="cell_class", values="mp")
        fig, ax = plt.subplots(figsize=(0.6 * m.shape[1] + 3.5, 0.32 * m.shape[0] + 2))
        im = ax.imshow(m.to_numpy(float), cmap="PuOr_r", vmin=0.3, vmax=0.7, aspect="auto")
        ax.set_xticks(range(m.shape[1]), m.columns, rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(m.shape[0]), m.index, fontsize=7)
        tiers = h.set_index(["gene_set", "cell_class"]).tier
        for i, gsn in enumerate(m.index):
            for j, cls in enumerate(m.columns):
                if tiers.get((gsn, cls), "") == "replicated":
                    ax.text(j, i, "*", ha="center", va="center", fontsize=8)
        ax.set_title(f"{ds}: mean variability percentile (0.5 = matched genes; * replicated)", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        out.figure(fig, f"donor_variability_{ds}", f"Between-donor variability of gene sets, {ds}")
        plt.close(fig)


if __name__ == "__main__":
    main()
