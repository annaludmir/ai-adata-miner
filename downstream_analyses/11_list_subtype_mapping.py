#!/usr/bin/env python3
"""11 - Which fine clusters, and which neuron sub-types, each gene list sits in.

Question: 04/06 place a list in a broad cell class. Within a class, is the
list concentrated in particular sub-types (fine clusters), consistently in
both donor sets? And among neurons, does the list follow the main sub-type
axes: excitatory vs inhibitory, and deep- vs upper-layer identity?

Method
  A. Cluster concentration. Cluster pseudobulks as in 07/08 (finest
     clustering, >= 100 cells, not donor-dominated; log2 TMM-CPM). Each gene
     is Z-scored across clusters; a set's score in a cluster is the mean Z
     of its genes. Null: random sets drawing each member from the same bin of
     mean level x spread (07's bins). Per cluster: effect in null SDs and an
     empirical p. v2 and v3 are combined per cluster (signed Stouffer, BH
     per dataset over set x cluster tests, tiered) -- which needs the same
     cluster in both chemistries. human_dev's cluster ids are shared (a
     label's v2 expression profile best matches the same label in v3; see
     _common.cluster_label_identity); cortex's fine clusterings
     were computed per chemistry (labels reused for different clusters), so
     there clusters are paired by expression: mutual best Pearson match of
     their gene-centred profiles over the 2,000 most variable genes, r >=
     0.5. Unpaired clusters are tested per chemistry only.
     Per set, the cell classes and regions of its replicated clusters are
     compared with those of all tested clusters.
  B. Neuron sub-type axes. Neuron clusters only (dominant class Neuron, purity
     >= 0.6). Axis score per cluster = mean Z of one marker group minus the
     other: excitatory (SLC17A6, SLC17A7, NEUROD2, NEUROD6) vs inhibitory
     (GAD1, GAD2, SLC32A1, DLX5), in human_dev (cortex is the EMX1 lineage,
     so has almost no inhibitory neurons); deep layer (BCL11B, TBR1, FEZF2,
     SOX5) vs upper layer (SATB2, CUX2, POU3F2, POU3F3), in cortex and in
     human_dev's telencephalic excitatory clusters. Test: Spearman of a set's
     score with the axis across neuron clusters, against matched random sets;
     the axis markers are removed from sets and from the null pool. v2 x v3
     combined.

Inputs (csv_exports/):
  <ds>__<chem>/09_pseudobulk/<finest clustering>__{pseudobulk_counts,group_summary}.csv
  <ds>__<chem>/04_clusters/cluster_profile_<clustering>.csv
  <ds>__<chem>/10_markers/top_markers_<finest clustering>.csv (labels only)
  <ds>__v2/11_panels/panel_coverage.csv; gene lists (config.GENE_LISTS_DIR)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "11_list_subtype_mapping"
TITLE = "Gene lists at sub-type resolution: fine clusters and neuron sub-type axes"
N_RANDOM = 1000
MIN_GENES = 5
MAX_LISTED = 6
NEURON = "Neuron"
MIN_PURITY = 0.6
MIN_AXIS_CLUSTERS = 12
EXC = ["SLC17A6", "SLC17A7", "NEUROD2", "NEUROD6"]
INH = ["GAD1", "GAD2", "SLC32A1", "DLX5"]
DEEP = ["BCL11B", "TBR1", "FEZF2", "SOX5"]
UPPER = ["SATB2", "CUX2", "POU3F2", "POU3F3"]
# axis name, (label A, genes A), (label B, genes B), datasets it applies to
AXES = [("excitatory_vs_inhibitory", ("excitatory", EXC), ("inhibitory", INH), ("human_dev",)),
        ("deep_vs_upper_layer", ("deep layer", DEEP), ("upper layer", UPPER), ("cortex", "human_dev"))]
TELENCEPHALIC = {"Telencephalon"}
MATCH_GENES = 2000
MATCH_MIN_R = 0.5


def gene_sets(out: C.Output) -> dict[str, dict[str, list[str]]]:
    return {ds: C.analysis_gene_sets(ds, groups=("ndd",), modules=False, out=out) for ds in C.DATASETS}


# ---------------------------------------------------------------------------
# A. cluster concentration
# ---------------------------------------------------------------------------
def cluster_scores(out: C.Output, sets_by_ds: dict, rng: np.random.Generator):
    rows, annots, exprs = [], {}, {}
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        grouping = C.COEXPR_GROUPING[ds]
        clustering = grouping.replace("cluster_", "", 1)
        out.used(f"{n}/09_pseudobulk/{grouping}__pseudobulk_counts.csv",
                 f"{n}/04_clusters/cluster_profile_{clustering}.csv")
        lc, annot = C.cluster_expression(n)
        annots[(ds, chem)], exprs[(ds, chem)] = annot, lc
        Z = C.zscore_rows(lc.to_numpy(float))
        bins = C.expression_bins(lc)
        pos = {g: i for i, g in enumerate(lc.index)}
        C.log(f"  {n}: {lc.shape[1]} clusters, {lc.shape[0]:,} genes")
        for name, genes in sets_by_ds[ds].items():
            idx = np.array([pos[g] for g in genes if g in pos])
            if idx.size < MIN_GENES:
                continue
            obs = Z[idx].mean(axis=0)
            null = C.set_mean_rows(Z, C.matched_sets(bins, idx, N_RANDOM, rng))
            mu, sd = null.mean(axis=0), null.std(axis=0)
            p = ((np.abs(null - mu) >= np.abs(obs - mu) - 1e-12).sum(axis=0) + 1) / (N_RANDOM + 1)
            rows.append(pd.DataFrame({
                "dataset": ds, "chemistry": chem, "gene_set": name, "cluster": lc.columns,
                "n_genes": idx.size, "score": obs,
                "effect_vs_null_sd": np.divide(obs - mu, sd, out=np.full_like(obs, np.nan), where=sd > 0),
                "perm_p": p}))
    per = pd.concat(rows, ignore_index=True)
    return per, annots, exprs


def pair_clusters(annots: dict, exprs: dict) -> tuple[dict, pd.DataFrame]:
    """{(ds, chem): {cluster: pair id}} and a table of how clusters were paired."""
    pairing, rows = {}, []
    for ds in C.DATASETS:
        a2, a3 = annots[(ds, "v2")], annots[(ds, "v3")]
        l2, l3 = exprs[(ds, "v2")], exprs[(ds, "v3")]
        shared = [c for c in l2.columns.intersection(l3.columns)]
        agree, _ = C.cluster_label_identity(ds, C.COEXPR_GROUPING[ds].replace("cluster_", "", 1))
        if agree >= C.SHARED_LABEL_IDENTITY:
            pairing[(ds, "v2")] = {c: c for c in shared}
            pairing[(ds, "v3")] = {c: c for c in shared}
            rows += [{"dataset": ds, "pair": c, "cluster_v2": c, "cluster_v3": c, "method": "shared label",
                      "r": np.nan} for c in shared]
            C.log(f"  {ds}: cluster labels shared across chemistries ({agree:.0%} best-match themselves)")
            continue
        genes = l2.index.intersection(l3.index)
        z2 = l2.loc[genes].sub(l2.loc[genes].mean(axis=1), axis=0)
        z3 = l3.loc[genes].sub(l3.loc[genes].mean(axis=1), axis=0)
        top = (z2.std(axis=1) + z3.std(axis=1)).sort_values(ascending=False).index[:MATCH_GENES]
        R = C.corr_rows(z2.loc[top].T.to_numpy(float), z3.loc[top].T.to_numpy(float))
        best3, best2 = R.argmax(axis=1), R.argmax(axis=0)
        p2, p3 = {}, {}
        for i, j in enumerate(best3):
            if best2[j] == i and R[i, j] >= MATCH_MIN_R:
                c2, c3 = l2.columns[i], l3.columns[j]
                pid = f"{c2}~{c3}"
                p2[c2], p3[c3] = pid, pid
                rows.append({"dataset": ds, "pair": pid, "cluster_v2": c2, "cluster_v3": c3,
                             "method": "expression match", "r": float(R[i, j]),
                             "class_v2": a2.loc[c2, "dominant_cell_class"],
                             "class_v3": a3.loc[c3, "dominant_cell_class"]})
        pairing[(ds, "v2")], pairing[(ds, "v3")] = p2, p3
        C.log(f"  {ds}: labels not shared ({agree:.0%} best-match themselves); {len(p2)} clusters paired "
              f"by expression (of {l2.shape[1]} v2 / {l3.shape[1]} v3)")
    return pairing, pd.DataFrame(rows)


def annotate(df: pd.DataFrame, annots: dict) -> pd.DataFrame:
    parts = []
    for ds, g in df.groupby("dataset", sort=False):
        mk = C.top_markers(C.ns(ds, "v2"), C.COEXPR_GROUPING[ds])
        for chem in C.CHEMISTRIES:
            a = annots[(ds, chem)].add_suffix(f"_{chem}")
            g = g.join(a, on=f"cluster_{chem}")
        g["dominant_cell_class"] = g["dominant_cell_class_v2"].fillna(g["dominant_cell_class_v3"])
        g["dominant_region"] = g["dominant_region_v2"].fillna(g["dominant_region_v3"])
        g["top_markers"] = g["cluster_v2"].astype(str).map(mk).fillna("")
        parts.append(g)
    return pd.concat(parts, ignore_index=True)


def concentration(comb: pd.DataFrame) -> pd.DataFrame:
    """Per set: classes / regions over-represented among its replicated enriched clusters."""
    rows = []
    for (ds, name), g in comb.groupby(["dataset", "gene_set"], sort=False):
        hit = g[(g.tier == "replicated") & (g.direction == "enriched")]
        for kind in ("dominant_cell_class", "dominant_region"):
            tested = g[kind].value_counts()
            hits = hit[kind].value_counts()
            for lvl, nt in tested.items():
                k = int(hits.get(lvl, 0))
                rows.append({"dataset": ds, "gene_set": name, "annotation": kind.replace("dominant_", ""),
                             "level": lvl, "n_clusters_tested": int(nt), "n_replicated_enriched": k,
                             "share_of_hits": k / len(hit) if len(hit) else np.nan,
                             "share_of_tested": nt / len(g),
                             "ratio": (k / len(hit)) / (nt / len(g)) if len(hit) else np.nan})
    return pd.DataFrame(rows)


def profile_agreement(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (ds, name), g in per.groupby(["dataset", "gene_set"], sort=False):
        w = g.dropna(subset=["pair"]).pivot_table(index="pair", columns="chemistry", values="score").dropna()
        if set(C.CHEMISTRIES) <= set(w.columns) and len(w) >= 10:
            rows.append({"dataset": ds, "gene_set": name, "n_shared_clusters": len(w),
                         "spearman_v2_v3": float(w["v2"].corr(w["v3"], method="spearman"))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# B. neuron sub-type axes
# ---------------------------------------------------------------------------
def axis_tests(sets_by_ds: dict, annots: dict, exprs: dict,
               rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, info = [], []
    for ds, chem in C.STRATA:
        lc, annot = exprs[(ds, chem)], annots[(ds, chem)]
        a = annot.loc[lc.columns]
        neuron = (a.dominant_cell_class.astype(str) == NEURON) & (a.purity_cell_class >= MIN_PURITY)
        for axis, (la, ga), (lb, gb), datasets in AXES:
            if ds not in datasets:
                continue
            cols = neuron.copy()
            if axis == "deep_vs_upper_layer" and ds == "human_dev":
                exc = lc.reindex([g for g in EXC if g in lc.index]).mean(axis=0)
                inh = lc.reindex([g for g in INH if g in lc.index]).mean(axis=0)
                cols &= a.dominant_region.astype(str).isin(TELENCEPHALIC) & (exc > inh)
            sub = lc.loc[:, cols.to_numpy()]
            A = [g for g in ga if g in sub.index]
            B = [g for g in gb if g in sub.index]
            row = {"dataset": ds, "chemistry": chem, "axis": axis, "n_clusters": sub.shape[1],
                   "markers_a": "|".join(A), "markers_b": "|".join(B)}
            if sub.shape[1] < MIN_AXIS_CLUSTERS or len(A) < 2 or len(B) < 2:
                info.append({**row, "status": "skipped (too few clusters or markers)"})
                continue
            Z = C.zscore_rows(sub.to_numpy(float))
            zi = {g: i for i, g in enumerate(sub.index)}
            sa, sb = Z[[zi[g] for g in A]].mean(axis=0), Z[[zi[g] for g in B]].mean(axis=0)
            ax = sa - sb
            info.append({**row, "status": "tested",
                         "spearman_a_vs_b_markers": float(pd.Series(sa).corr(pd.Series(sb), method="spearman"))})
            keep = np.array([g not in set(A) | set(B) for g in sub.index])
            Zk, genes_k = Z[keep], sub.index[keep]
            bins = C.expression_bins(sub.loc[keep])
            pos = {g: i for i, g in enumerate(genes_k)}
            for name, genes in sets_by_ds[ds].items():
                idx = np.array([pos[g] for g in genes if g in pos])
                if idx.size < MIN_GENES:
                    continue
                rho = float(C.spearman_rows(Zk[idx].mean(axis=0)[None, :], ax)[0])
                null = C.spearman_rows(C.set_mean_rows(Zk, C.matched_sets(bins, idx, N_RANDOM, rng)), ax)
                eff, p, mu, sd = C.null_effect(rho, null)
                rows.append({"dataset": ds, "chemistry": chem, "axis": axis, "gene_set": name,
                             "pole_a": la, "pole_b": lb, "n_genes": int(idx.size),
                             "n_clusters": sub.shape[1], "rho_with_axis": rho, "null_mean": mu,
                             "null_sd": sd, "effect_vs_null_sd": eff, "perm_p": p})
    return pd.DataFrame(rows), pd.DataFrame(info)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sets_by_ds = gene_sets(out)

    per, annots, exprs = cluster_scores(out, sets_by_ds, rng)
    pairing, pairs = pair_clusters(annots, exprs)
    per["pair"] = [pairing[(d, c)].get(cl) for d, c, cl in zip(per.dataset, per.chemistry, per.cluster)]
    out.write(pairs, "cluster_pairs",
              "How clusters were matched between v2 and v3: shared label, or mutual best "
              "expression match (r = Pearson of gene-centred profiles)")
    for ds, chem in C.STRATA:
        if (C.EXPORTS / C.ns(ds, chem) / "10_markers" / f"top_markers_{C.COEXPR_GROUPING[ds]}.csv").exists():
            out.used(f"{C.ns(ds, chem)}/10_markers/top_markers_{C.COEXPR_GROUPING[ds]}.csv")
    comb = C.combine_chemistries(per.dropna(subset=["pair"]), ["dataset", "gene_set", "pair"],
                                 labels=("enriched", "depleted"), carry=("score", "cluster"))
    comb = annotate(comb, annots)
    out.write(per, "cluster_scores_per_stratum",
              "Per set x cluster x stratum: mean Z of set genes, effect vs matched random sets, p")
    keep_cols = ["dataset", "gene_set", "pair", "cluster_v2", "cluster_v3", "dominant_cell_class", "dominant_region",
                 "top_markers", "score_v2", "score_v3", "effect_v2", "p_v2", "effect_v3", "p_v3",
                 "stouffer_z", "combined_p", "combined_q", "tier", "direction"]
    keep_cols += [c for c in comb.columns if c.startswith(("n_cells_", "age_pcw_median_", "purity_cell_class_"))]
    out.write(comb[keep_cols], "cluster_scores_combined",
              "Per set x cluster: v2/v3 combined enrichment, tier, cluster annotation and top markers")
    conc = concentration(comb)
    out.write(conc, "list_cluster_concentration",
              "Per set: share of its replicated enriched clusters in each class / region vs share "
              "of tested clusters (ratio > 1 = over-represented)")
    agree = profile_agreement(per)
    out.write(agree, "cluster_profile_agreement",
              "Per set: Spearman of its cluster scores between v2 and v3 over paired clusters")

    axes, axis_info = axis_tests(sets_by_ds, annots, exprs, rng)
    out.write(axis_info, "neuron_axes", "Neuron sub-type axes: clusters and markers used per stratum")
    acomb = pd.DataFrame()
    if not axes.empty:
        out.write(axes, "neuron_axis_tests_per_stratum",
                  "Per set x axis x stratum: Spearman of set score with the axis across neuron "
                  "clusters, vs matched random sets")
        acomb = C.combine_chemistries(axes, ["dataset", "axis", "gene_set"],
                                      labels=("pole a", "pole b"), carry=("rho_with_axis", "n_genes"))
        poles = axes.drop_duplicates(["axis"]).set_index("axis")[["pole_a", "pole_b"]]
        acomb["direction"] = [poles.loc[r.axis, "pole_a"] if r.stouffer_z > 0 else poles.loc[r.axis, "pole_b"]
                              for r in acomb.itertuples()]
        out.write(acomb, "neuron_axis_tests_combined",
                  "Per set x axis: v2/v3 combined association with the axis; direction = the pole "
                  "the set leans to")

    figures(out, conc, acomb)

    # ---- findings -----------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        c = comb[comb.dataset == ds]
        n_cl = c.pair.nunique()
        hits = c[(c.tier == "replicated") & (c.direction == "enriched")]
        lists = [s for s in c.gene_set.unique() if s.startswith("list:")] or list(c.gene_set.unique())
        parts = []
        for name in lists:
            h = hits[hits.gene_set == name].sort_values("stouffer_z", ascending=False)
            cc = conc[(conc.dataset == ds) & (conc.gene_set == name) & (conc.annotation == "cell_class")
                      & (conc.n_replicated_enriched >= 2) & (conc.ratio >= 1.5)].sort_values("ratio", ascending=False)
            if h.empty:
                parts.append(f"{name}: no replicated cluster")
                continue
            top = "; ".join(f"cluster {r.pair} ({r.dominant_cell_class}"
                            + (f", {r.dominant_region}" if ds == "human_dev" else "")
                            + (f"; {r.top_markers.split(', ')[0]}, {r.top_markers.split(', ')[1]}"
                               if r.top_markers.count(",") >= 1 else "")
                            + f"; {r.effect_v2:+.1f}/{r.effect_v3:+.1f})" for r in h.head(3).itertuples())
            over = (" -- over-represented: " + ", ".join(f"{r.level} ({r.ratio:.1f}x)" for r in cc.head(3).itertuples())
                    if len(cc) else "")
            parts.append(f"{name}: {len(h)} replicated clusters, top {top}{over}")
        how = (pairs[pairs.dataset == ds].method.iloc[0] if (pairs.dataset == ds).any() else "none")
        f.append(f"**{ds}: clusters where each set concentrates** ({n_cl} clusters in both "
                 f"chemistries, paired by {how}; effect in null SDs v2/v3): " + " | ".join(parts[:MAX_LISTED * 2]) + ".")
    if not agree.empty:
        f.append("**Cluster profiles reproduce across donor sets** (Spearman of set scores v2 vs v3 "
                 "over paired clusters): " + "; ".join(
                     f"{ds} median {g.spearman_v2_v3.median():.2f} (range {g.spearman_v2_v3.min():.2f}"
                     f"-{g.spearman_v2_v3.max():.2f})" for ds, g in agree.groupby("dataset")) + ".")
    if not acomb.empty:
        for (ds, axis), g in acomb.groupby(["dataset", "axis"], sort=False):
            inf = axis_info[(axis_info.dataset == ds) & (axis_info.axis == axis)]
            ncl = "/".join(str(int(x)) for x in inf.n_clusters)
            rep = g[g.tier != ""]
            f.append(f"**{ds}, {axis.replace('_', ' ')}** ({ncl} neuron clusters v2/v3; marker groups "
                     f"correlate rho {inf.spearman_a_vs_b_markers.round(2).tolist()}): "
                     + ("; ".join(f"{r.gene_set} toward {r.direction} (rho {r.rho_with_axis_v2:+.2f}/"
                                  f"{r.rho_with_axis_v3:+.2f}; {r.tier})" for r in rep.head(MAX_LISTED * 2).itertuples())
                        if len(rep) else "no set leans to either pole in both donor sets") + ".")
    skipped = axis_info[axis_info.status != "tested"] if not axis_info.empty else axis_info
    if len(skipped):
        f.append("**Axes not tested** (too few neuron clusters or markers): "
                 + "; ".join(f"{r.dataset} {r.chemistry} {r.axis} ({r.n_clusters} clusters)"
                             for r in skipped.itertuples()) + ".")

    out.summary(
        TITLE,
        "Within broad cell classes, which fine clusters (sub-types) does each gene list concentrate "
        "in, consistently in both donor sets; and among neurons, does the list follow the "
        "excitatory-inhibitory and deep-upper layer axes?",
        ["A: cluster pseudobulks (07/08 rules), gene Z-scores across clusters, set score = mean Z; "
         f"{N_RANDOM:,} random sets matched member by member on mean level x spread; per cluster "
         "effect in null SDs and empirical p; the same cluster in v2 and v3 (shared label, or "
         f"mutual best expression match with r >= {MATCH_MIN_R} where labels differ between "
         "chemistries) combined by signed Stouffer, BH per dataset over set x cluster tests, tiered. Classes / regions of a set's "
         "replicated clusters compared with all tested clusters.",
         f"B: neuron clusters (purity >= {MIN_PURITY}); axis = mean Z of one marker group minus "
         f"the other (excitatory {', '.join(EXC)} vs inhibitory {', '.join(INH)}, human_dev; deep "
         f"{', '.join(DEEP)} vs upper {', '.join(UPPER)}, cortex and human_dev telencephalic "
         "excitatory clusters). Spearman of set score with the axis vs matched random sets, axis "
         "markers removed; v2 x v3 combined.",
         "Seed NDD panels and user lists (GWAS lists one gene per locus)."],
        f,
        ["Clusters are not donors: a cluster's v2 and v3 cells come from different people, which "
         "is what makes the per-cluster replication meaningful, but within a chemistry a cluster "
         "can be dominated by few donors (one-donor clusters are excluded).",
         "Cluster scores are relative to the other clusters of the same file; a list 'enriched' in "
         "a cluster is high there compared with the rest of the dataset, not necessarily specific.",
         "Neuron axes are confounded with maturation: deep-layer neurons are born first, so a "
         "list leaning deep may simply be higher in older neurons."],
        ["Pair with B5 (pseudotime) to separate sub-type identity from maturation.",
         "Name clusters with a reference atlas instead of their top markers."])


def figures(out: C.Output, conc: pd.DataFrame, acomb: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None:
        return
    for ds in C.DATASETS:
        h = conc[(conc.dataset == ds) & (conc.annotation == "cell_class")]
        if h.empty:
            continue
        frac = h.assign(f=h.n_replicated_enriched / h.n_clusters_tested).pivot(
            index="gene_set", columns="level", values="f")
        fig, ax = plt.subplots(figsize=(0.55 * frac.shape[1] + 3.5, 0.35 * frac.shape[0] + 2))
        im = ax.imshow(frac.to_numpy(float), cmap="Reds", vmin=0, vmax=max(0.3, np.nanmax(frac.to_numpy(float))),
                       aspect="auto")
        ax.set_xticks(range(frac.shape[1]), frac.columns, rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(frac.shape[0]), frac.index, fontsize=7)
        ax.set_title(f"{ds}: share of each class's clusters where the set is enriched (replicated)", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        out.figure(fig, f"cluster_concentration_{ds}", f"Replicated cluster enrichment by cell class, {ds}")
        plt.close(fig)
    if acomb.empty:
        return
    groups = list(acomb.groupby(["dataset", "axis"], sort=False))
    fig, axes = plt.subplots(1, len(groups), figsize=(4.2 * len(groups), 0.3 * acomb.gene_set.nunique() + 2),
                             squeeze=False)
    for ax, ((ds, axis), g) in zip(axes.flat, groups):
        g = g.sort_values("stouffer_z")
        colors = ["#c05621" if t else "#a0aec0" for t in g.tier]
        ax.barh(range(len(g)), g.stouffer_z, color=colors)
        ax.set_yticks(range(len(g)), g.gene_set, fontsize=6)
        ax.axvline(0, color="grey", lw=0.6)
        ax.set_title(f"{ds}: {axis.replace('_', ' ')}\n(+ = toward first pole; orange = tiered)", fontsize=8)
        ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    out.figure(fig, "neuron_axes", "Gene sets along neuron sub-type axes (combined Z)")
    plt.close(fig)


if __name__ == "__main__":
    main()
