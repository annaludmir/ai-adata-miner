#!/usr/bin/env python3
"""18 - Co-expression rewiring: do gene lists become more or less coherent with age?

Question: 07 measured whether a list's genes co-express across clusters. Does
that coherence change between early and late development -- a programme
assembling (rising coherence) or loosening -- consistently in both donor sets?

Method
  Cluster pseudobulks as in 07/08 (finest clustering, >= 100 cells, not
  donor-dominated; log2 TMM-CPM). Within each dominant cell class with >= 10
  clusters, clusters are split at the class's median cluster age (04's median
  age of the cluster's cells) into an early and a late half, so both halves
  hold the same mix of classes. Coherence = mean pairwise Spearman of the
  set's genes over the clusters of one half, in two contexts as in 07:
  across clusters, and within class (each gene's class mean removed inside
  the half). Statistic: late minus early coherence. Null: random sets drawn
  member by member from the same bin of mean level x spread (07's bins, from
  all clusters). v2 x v3: signed Stouffer, BH per dataset x context, tiered.

Inputs (csv_exports/):
  <ds>__<chem>/09_pseudobulk/<finest clustering>__{pseudobulk_counts,group_summary}.csv
  <ds>__<chem>/04_clusters/cluster_profile_<clustering>.csv
  <ds>__v2/11_panels/panel_coverage.csv; gene lists; results/08 modules
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "18_coexpression_rewiring"
TITLE = "Co-expression rewiring: coherence of gene lists in early vs late clusters"
N_RANDOM = 500
MIN_GENES = 5
MIN_CLASS_CLUSTERS = 10
MAX_LISTED = 8


def split_by_age(annot: pd.DataFrame, cols: list[str]) -> tuple[list[str], list[str], dict]:
    a = annot.loc[cols]
    early, late, info = [], [], {}
    for cls, g in a.groupby("dominant_cell_class"):
        if len(g) < MIN_CLASS_CLUSTERS:
            continue
        order = g.sort_values(["age_pcw_median", "n_cells"]).index.tolist()
        half = len(order) // 2
        early += order[:half]
        late += order[len(order) - half:]
        info[cls] = (half, float(g.loc[order[:half], "age_pcw_median"].median()),
                     float(g.loc[order[len(order) - half:], "age_pcw_median"].median()))
    return early, late, info


def half_values(lc: pd.DataFrame, annot: pd.DataFrame, cols: list[str], context: str) -> np.ndarray:
    X = lc[cols].to_numpy(float)
    if context == "within_class":
        cls = annot.loc[cols, "dominant_cell_class"].astype(str).to_numpy()
        X = X.copy()
        for c in np.unique(cls):
            m = cls == c
            X[:, m] -= X[:, m].mean(axis=1, keepdims=True)
    return C.unit_rank_rows(X)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sets = {ds: C.analysis_gene_sets(ds, out=out) for ds in C.DATASETS}
    rows, splits = [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        grouping = C.COEXPR_GROUPING[ds]
        out.used(f"{n}/09_pseudobulk/{grouping}__pseudobulk_counts.csv",
                 f"{n}/04_clusters/cluster_profile_{grouping.replace('cluster_', '', 1)}.csv")
        lc, annot = C.cluster_expression(n)
        early, late, info = split_by_age(annot, list(lc.columns))
        for cls, (k, ae, al) in info.items():
            splits.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "clusters_per_half": k,
                           "median_age_early": ae, "median_age_late": al})
        if len(early) < 2 * MIN_CLASS_CLUSTERS // 2:
            continue
        bins = C.expression_bins(lc)
        pos = {g: i for i, g in enumerate(lc.index)}
        for context in ("across_clusters", "within_class"):
            Ue = half_values(lc, annot, early, context)
            Ul = half_values(lc, annot, late, context)
            for name, genes in sets[ds].items():
                idx = np.array([pos[g] for g in genes if g in pos])
                if idx.size < MIN_GENES:
                    continue
                k = idx.size
                coh_e = float(C.coherence_from_sums(Ue[idx].sum(axis=0)[None, :], k)[0])
                coh_l = float(C.coherence_from_sums(Ul[idx].sum(axis=0)[None, :], k)[0])
                rand = C.matched_sets(bins, idx, N_RANDOM, rng)
                null = (C.coherence_from_sums(C.set_sums(Ul, rand), k)
                        - C.coherence_from_sums(C.set_sums(Ue, rand), k))
                eff, p, mu, sd = C.null_effect(coh_l - coh_e, null)
                rows.append({"dataset": ds, "chemistry": chem, "context": context, "gene_set": name,
                             "n_genes": k, "n_clusters_per_half": len(early), "coherence_early": coh_e,
                             "coherence_late": coh_l, "late_minus_early": coh_l - coh_e,
                             "null_mean": mu, "effect_vs_null_sd": eff, "perm_p": p})
    per = pd.DataFrame(rows)
    splits = pd.DataFrame(splits)
    comb = pd.concat([C.combine_chemistries(g, ["dataset", "context", "gene_set"],
                                            labels=("more coherent late", "less coherent late"),
                                            carry=("coherence_early", "coherence_late"))
                      for _, g in per.groupby("context")], ignore_index=True) if not per.empty else per
    out.write(splits, "age_split", "Per stratum x class: clusters per half and median cluster age of each half")
    out.write(per, "rewiring_per_stratum", "Per stratum x context x set: coherence in early and late "
              "clusters, late - early vs matched random sets")
    out.write(comb, "rewiring_combined", "v2 x v3 combined change in coherence; BH per dataset within "
              "context; tier")
    figures(out, comb)

    # ---- findings -----------------------------------------------------------
    f = []
    if not splits.empty:
        f.append("**Age split** (median cluster age early -> late, per class): "
                 + "; ".join(f"{ds} {chem}: " + ", ".join(
                     f"{r.cell_class} {r.median_age_early:g}->{r.median_age_late:g}"
                     for r in g.itertuples()) for (ds, chem), g in splits.groupby(["dataset", "chemistry"]))
                 + ".")
    for ds in C.DATASETS:
        for context in ("across_clusters", "within_class"):
            g = comb[(comb.dataset == ds) & (comb.context == context)] if not comb.empty else comb
            rep = g[g.tier != ""] if len(g) else g
            f.append(f"**{ds}, {context.replace('_', ' ')}: sets whose coherence changes with age** "
                     "(coherence early -> late, v2 / v3): "
                     + ("; ".join(f"{r.gene_set} {r.direction} ({r.coherence_early_v2:.2f}->{r.coherence_late_v2:.2f}"
                                  f" / {r.coherence_early_v3:.2f}->{r.coherence_late_v3:.2f}; {r.tier})"
                                  for r in rep.head(MAX_LISTED).itertuples())
                        + ("" if len(rep) <= MAX_LISTED else f" (+{len(rep) - MAX_LISTED} more)")
                        if len(rep) else "none replicated") + ".")
    out.summary(
        TITLE,
        "Do gene lists and modules co-express more or less tightly in late than in early "
        "clusters, within the same mix of cell classes, in both donor sets?",
        [f"Cluster pseudobulks (07/08 rules). Per class with >= {MIN_CLASS_CLUSTERS} clusters, clusters "
         "split at the class's median cluster age; halves pooled over classes.",
         "Coherence = mean pairwise Spearman over one half's clusters; contexts across clusters and "
         f"within class (class means removed inside each half). Late - early vs {N_RANDOM} random sets "
         "matched on mean level x spread; v2 x v3 signed Stouffer, BH per dataset x context, tiered."],
        f,
        ["A cluster's age is the median age of its cells; clusters pooling several donors blur the "
         "split, which biases changes towards zero.",
         "Early and late halves can differ in sub-type mix within a class, which changes "
         "co-expression without any rewiring inside cells.",
         "Coherence over ~half the clusters is noisier than 07's estimate over all of them."],
        ["B4 (within-cell co-expression per age) would test rewiring inside one cell type directly."])


def figures(out: C.Output, comb: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or comb.empty:
        return
    h = comb[comb.gene_set.str.startswith(("list:", "seed:"))]
    if h.empty:
        return
    groups = list(h.groupby(["dataset", "context"], sort=False))
    fig, axes = plt.subplots(1, len(groups), figsize=(3.6 * len(groups), 0.3 * h.gene_set.nunique() + 2),
                             squeeze=False)
    for ax, ((ds, ctx), g) in zip(axes.flat, groups):
        g = g.sort_values("stouffer_z")
        ax.barh(range(len(g)), g.stouffer_z, color=["#c05621" if t else "#a0aec0" for t in g.tier])
        ax.set_yticks(range(len(g)), g.gene_set, fontsize=6)
        ax.axvline(0, color="grey", lw=0.6)
        ax.set_title(f"{ds}, {ctx.replace('_', ' ')}\n(+ = more coherent late)", fontsize=8)
        ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    out.figure(fig, "coexpression_rewiring", "Change in gene-set coherence, late vs early clusters (combined Z)")
    plt.close(fig)


if __name__ == "__main__":
    main()
