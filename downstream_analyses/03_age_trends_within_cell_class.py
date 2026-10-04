#!/usr/bin/env python3
"""03 - Which genes change with age inside a cell class, in both donor sets?

Question: holding cell identity fixed, which genes rise or fall across
development -- and do those trends replicate between the two chemistries,
which are disjoint sets of donors?

Working inside a cell class matters: a gene can look "up with age" across all
cells purely because neurons, which express it, become more common (02 shows
they do). Within-class means remove that composition effect.

Method: for each stratum x cell class, pseudobulk counts per age point are
TMM-normalised to log2 CPM and every expressed gene is correlated with age by
Spearman; with 5-9 age
points the p-value comes from an exact permutation over all age orderings
(Monte Carlo above 9). Each age point is one donor (occasionally 2-3), so n is
the number of age points -- never cells. v2 and v3 are combined by signed
Stouffer and tiered as in 02.

Why TMM log CPM and not mean_lognorm: in cortex, UMIs per cell fall with age
in every class, and mean log1p(CP10K) inherits that as a gene-wide downward
drift (median gene rho ~ -0.5, i.e. most genes "decrease"). TMM removes it;
the residual median rho per class is written to trend_counts.csv as a check.

Primary result: cortex. human_dev classes pool every brain region, and which
regions were dissected changes with age, so its within-class trends mix
development with regional sampling; they are reported, flagged.

Inputs (csv_exports/):
  <ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,detection_fraction,group_summary}.csv
  <ds>__<chem>/11_panels/panel_coverage.csv
  _cross_dataset/gene_id_map.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "03_age_trends_within_cell_class"
MIN_AGES = 5            # age points per stratum x class
MIN_DETECTION = 0.10    # gene must be detected in >= 10% of cells ...
MIN_DETECTED_AGES = 2   # ... in at least this many age points
MIN_MEAN_CPM = 5.0      # and average >= 5 CPM across age points
N_PERM = 400_000        # exact up to 9 ages (9! = 362,880)
TOP_N = 15              # genes listed per class and direction in the summary
PRIMARY = "cortex"
DRIFT_RHO = 0.1        # |median gene rho| above this = a global shift remains


def test_stratum(n: str, ds: str, chem: str, chrom, panel_of) -> pd.DataFrame:
    rows = []
    for cls, (ln, det, ages, min_cells) in C.class_age_logcpm(n).items():
        if len(ages) < MIN_AGES:
            C.log(f"  {n} {cls}: {len(ages)} age points (<{MIN_AGES}), skipped")
            continue
        expressed = (((det >= MIN_DETECTION).sum(axis=1) >= MIN_DETECTED_AGES)
                     & (np.log2(MIN_MEAN_CPM + 1) <= ln.mean(axis=1)))
        X = ln.loc[expressed].to_numpy(float)
        genes = ln.index[expressed].astype(str)
        rho = C.spearman_rows(X, ages)
        p, exact = C.spearman_perm_p(rho, ages, n_perm=N_PERM)
        xc = ages - ages.mean()
        slope = (X - X.mean(axis=1, keepdims=True)) @ xc / (xc ** 2).sum()
        C.log(f"  {n} {cls}: {len(genes):,} expressed genes x {len(ages)} ages "
              f"({ages.min():g}-{ages.max():g} pcw), {'exact' if exact else 'Monte Carlo'} null")
        rows.append(pd.DataFrame({
            "dataset": ds, "chemistry": chem, "cell_class": cls, "gene": genes,
            "n_ages": len(ages), "age_range": f"{ages.min():g}-{ages.max():g}",
            "min_cells_per_age": min_cells,
            "spearman_rho": rho, "perm_p": p, "exact": exact,
            "log2fc_per_week": slope,
            "mean_log2cpm": X.mean(axis=1), "max_detection": det.loc[expressed].max(axis=1).to_numpy(),
            "sex_linked": C.sex_linked(genes, chrom),
            "panels": [panel_of.get(g, "") for g in genes],
        }))
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not df.empty:
        df["q_within_stratum"] = np.nan
        for _, idx in df.groupby("cell_class").groups.items():
            df.loc[idx, "q_within_stratum"] = C.bh(df.loc[idx, "perm_p"])
    return df


def strongest(g: pd.DataFrame, direction: str, n: int) -> list[str]:
    """Top-n genes in one direction, ranked reproducibly.

    Exact permutation p-values are discrete, so many genes share a Stouffer Z
    and the order among them would otherwise depend on floating-point noise
    (it differed between the cluster and a laptop at the 1e-15 level). Ties
    are broken by effect size -- mean |log2 FC per week| over both
    chemistries -- then by gene name.
    """
    h = g[g.direction == direction].assign(
        _z=lambda d: d.stouffer_z.abs().round(9),
        _eff=lambda d: (d.log2fc_per_week_v2.abs() + d.log2fc_per_week_v3.abs()) / 2)
    return list(h.sort_values(["_z", "_eff", "gene"], ascending=[False, False, True])
                .gene.head(n))


def combine(per: pd.DataFrame) -> pd.DataFrame:
    """v2 x v3 per (dataset, class, gene) -> signed Stouffer, BH, tier."""
    out = []
    for (ds, cls), g in per.groupby(["dataset", "cell_class"], sort=False):
        v2 = g[g.chemistry == "v2"].set_index("gene")
        v3 = g[g.chemistry == "v3"].set_index("gene")
        both = v2.index.intersection(v3.index)
        if both.empty:
            continue
        a, b = v2.loc[both], v3.loc[both]
        Z, pc = C.signed_stouffer([a.spearman_rho.to_numpy(), b.spearman_rho.to_numpy()],
                                  [a.perm_p.to_numpy(), b.perm_p.to_numpy()],
                                  [np.sqrt(a.n_ages.iloc[0]), np.sqrt(b.n_ages.iloc[0])])
        q = C.bh(pc)
        df = pd.DataFrame({
            "dataset": ds, "cell_class": cls, "gene": both,
            "rho_v2": a.spearman_rho.to_numpy(), "p_v2": a.perm_p.to_numpy(),
            "n_ages_v2": a.n_ages.to_numpy(), "ages_v2": a.age_range.to_numpy(),
            "rho_v3": b.spearman_rho.to_numpy(), "p_v3": b.perm_p.to_numpy(),
            "n_ages_v3": b.n_ages.to_numpy(), "ages_v3": b.age_range.to_numpy(),
            "log2fc_per_week_v2": a.log2fc_per_week.to_numpy(),
            "log2fc_per_week_v3": b.log2fc_per_week.to_numpy(),
            "mean_log2cpm": (a.mean_log2cpm.to_numpy() + b.mean_log2cpm.to_numpy()) / 2,
            "stouffer_z": Z, "combined_p": pc, "combined_q": q,
            "sex_linked": a.sex_linked.to_numpy(), "panels": a.panels.to_numpy(),
        })
        df["tier"] = C.replication_tier(df.rho_v2, df.p_v2, df.rho_v3, df.p_v3, df.combined_q)
        df["direction"] = np.where(df.stouffer_z > 0, "up", "down")
        out.append(df)
    return pd.concat(out, ignore_index=True).sort_values(
        ["dataset", "cell_class", "combined_p"]) if out else pd.DataFrame()


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    out.used("_cross_dataset/gene_id_map.csv")
    chrom = C.chromosome_map()

    per = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv",
                 f"{n}/09_pseudobulk/cell_class_x_age__detection_fraction.csv",
                 f"{n}/09_pseudobulk/cell_class_x_age__group_summary.csv",
                 f"{n}/11_panels/panel_coverage.csv")
        pan = C.panels(n)
        pan = pan[pan.panel_group == "ndd"]
        panel_of = pan.groupby("gene")["panel"].apply(lambda s: "|".join(sorted(set(s)))).to_dict()
        per.append(test_stratum(n, ds, chem, chrom, panel_of))
    per = pd.concat(per, ignore_index=True)
    out.write(per, "age_trends_per_stratum",
              "Every expressed gene x class x stratum: Spearman vs age, exact permutation p")

    comb = combine(per)
    out.write(comb, "age_trends_combined",
              "v2 and v3 combined per gene x class; tier replicated / supported")

    rep = comb[comb.tier == "replicated"]
    counts = (comb.assign(n=1).pivot_table(index=["dataset", "cell_class"], columns="tier",
                                            values="n", aggfunc="sum", fill_value=0)
              .rename(columns={"": "not_replicated"}).reset_index())
    ups = rep[rep.direction == "up"].groupby(["dataset", "cell_class"]).size().rename("replicated_up")
    downs = rep[rep.direction == "down"].groupby(["dataset", "cell_class"]).size().rename("replicated_down")
    counts = counts.merge(ups.reset_index(), how="left").merge(downs.reset_index(), how="left").fillna(0)
    drift = (per.groupby(["dataset", "cell_class", "chemistry"])["spearman_rho"].median()
             .unstack("chemistry").add_prefix("median_gene_rho_").reset_index())
    counts = counts.merge(drift, how="left")
    out.write(counts, "trend_counts",
              "Genes per class by tier and direction; median_gene_rho_* near 0 = no global drift")

    ndd = comb[(comb.panels != "") & comb.tier.isin(["replicated", "supported"])]
    out.write(ndd, "ndd_genes_with_age_trends",
              "NDD-panel genes whose within-class age trend is replicated or supported")

    conc = concordance(comb)
    if not conc.empty:
        out.write(conc, "cortex_vs_human_dev_concordance",
                  "Correlation of combined Z between files, shared classes (same donors)")

    figures(out, per, rep)

    # ---- findings -----------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        c = counts[counts.dataset == ds]
        if c.empty:
            continue
        label = "primary" if ds == PRIMARY else "region-confounded, see limitations"
        f.append(f"**{ds} ({label})**: replicated within-class age trends per class -- "
                 + "; ".join(f"{r.cell_class} {int(r.replicated_up)} up / {int(r.replicated_down)} down"
                             f" of {int(r.not_replicated + r.get('replicated', 0) + r.get('supported', 0))} tested"
                             for _, r in c.iterrows()) + ".")
    for cls, g in rep[(rep.dataset == PRIMARY) & ~rep.sex_linked].groupby("cell_class"):
        dn = ", ".join(strongest(g, "down", TOP_N))
        up = ", ".join(strongest(g, "up", TOP_N))
        f.append(f"**{PRIMARY} {cls}** -- strongest replicated, rising: {up or 'none'}; "
                 f"falling: {dn or 'none'}.")
    drifted = counts[(counts.filter(like="median_gene_rho_").abs() > DRIFT_RHO).any(axis=1)]
    if len(drifted):
        f.append("**Residual global drift after TMM** (median gene rho beyond "
                 f"+/-{DRIFT_RHO}) in "
                 + "; ".join(f"{r.dataset} {r.cell_class} (v2 {r.median_gene_rho_v2:+.2f}, "
                             f"v3 {r.median_gene_rho_v3:+.2f})" for r in drifted.itertuples())
                 + ". Trends there are partly a whole-transcriptome shift (in erythrocytes, "
                 "haemoglobin taking over the transcriptome as they mature) -- prefer genes "
                 "whose |rho| clearly exceeds that offset.")
    sx = rep[rep.sex_linked]
    if len(sx):
        f.append(f"**{len(sx)} replicated trends involve sex-linked genes** "
                 f"({', '.join(sorted(set(sx.gene)))}); they track which donors sit at "
                 "which ages and are excluded from the lists above.")
    nd = ndd[(ndd.dataset == PRIMARY) & (ndd.tier == "replicated")]
    if len(nd):
        f.append(f"**NDD-panel genes with replicated trends in {PRIMARY}**: "
                 + "; ".join(f"{r.gene} {r.direction} in {r.cell_class} ({r.panels})"
                             for r in nd.sort_values("combined_p").head(25).itertuples())
                 + (" ..." if len(nd) > 25 else "") + ". Full list: ndd_genes_with_age_trends.csv.")
    if not conc.empty:
        f.append("**cortex and human_dev agree on within-class trends** (Spearman of combined Z, "
                 "shared classes): " + "; ".join(f"{r.cell_class} {r.spearman_z:.2f}"
                                                for r in conc.itertuples())
                 + ". Same donors, so this measures how much whole-brain pooling and processing "
                 "change the answer -- not independent replication.")

    out.summary(
        "Within-cell-class gene expression across age",
        "Holding cell identity fixed, which genes rise or fall across development, "
        "consistently in two independent donor sets?",
        ["Input: per-(cell class, age) pseudobulk counts from 09_pseudobulk, TMM-normalised "
         "within each class to log2 CPM (mean log1p(CP10K) drifts with UMIs per cell, which "
         f"falls with age in cortex). Age points with < {C.MIN_CELLS} cells dropped; classes "
         f"need >= {MIN_AGES} age points per chemistry.",
         f"Genes tested if detected in >= {MIN_DETECTION:.0%} of cells at >= {MIN_DETECTED_AGES} "
         f"age points and averaging >= {MIN_MEAN_CPM:g} CPM.",
         "Replicate unit: the age point (one donor, occasionally 2-3 at the same age). "
         "Spearman rho vs age; two-sided p exact over all age orderings up to 9 points, "
         f"Monte Carlo ({N_PERM:,}) above.",
         "v2 and v3 combined by signed Stouffer (weights sqrt(n ages)), BH per dataset x class. "
         "Replicated = same direction, q < 0.05, each chemistry nominal (one-sided p < 0.05).",
         "Sex-linked genes (chrY, XIST/TSIX) flagged and excluded from gene lists."],
        f,
        ["Each age point is effectively one donor: a trend is a trend across 5-9 people per "
         "chemistry. A gene specific to one unusual donor can still pass if that donor sits "
         "at an extreme age.",
         "v2 covers ~6.9-10 pcw and v3 ~5-14 pcw; replicated means monotonic across both "
         "windows, so transient (rise-then-fall) programmes are missed by design.",
         "Within-class means still mix sub-types: a 'Neuron' trend can be a shift in which "
         "neuron types are present (e.g. more deep-layer vs upper-layer neurons with age).",
         "human_dev classes pool regions whose dissection changes with age (see 01); treat "
         "its trends as hypotheses, and cortex as the cleaner test.",
         "Spearman captures monotonic trends; with 5-9 points, small effects are underpowered."],
        ["Repeat at cluster resolution within a class, to separate sub-type shifts from "
         "within-cell-type change.",
         "Add donors at ages each chemistry lacks; with >1 donor per age, model donor as random.",
         "Check top genes against an external developmental atlas (e.g. BrainSpan)."])


def concordance(comb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cls in sorted(set(comb[comb.dataset == "cortex"].cell_class)
                      & set(comb[comb.dataset == "human_dev"].cell_class)):
        a = comb[(comb.dataset == "cortex") & (comb.cell_class == cls)].set_index("gene")
        b = comb[(comb.dataset == "human_dev") & (comb.cell_class == cls)].set_index("gene")
        both = a.index.intersection(b.index)
        if len(both) < 50:
            continue
        za, zb = a.loc[both, "stouffer_z"], b.loc[both, "stouffer_z"]
        ra = a.loc[both, "tier"] == "replicated"
        rb = b.loc[both, "tier"] == "replicated"
        rows.append({"cell_class": cls, "n_shared_genes": len(both),
                     "spearman_z": za.corr(zb, method="spearman"),
                     "replicated_cortex": int(ra.sum()), "replicated_human_dev": int(rb.sum()),
                     "replicated_both_same_direction":
                         int((ra & rb & (np.sign(za) == np.sign(zb))).sum())})
    return pd.DataFrame(rows)


def figures(out: C.Output, per: pd.DataFrame, rep: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or rep.empty:
        return
    for cls, g in rep[(rep.dataset == PRIMARY) & ~rep.sex_linked].groupby("cell_class"):
        ups = strongest(g, "up", 4)
        downs = strongest(g, "down", 4)
        if not ups and not downs:
            continue
        slots = ups + [None] * (4 - len(ups)) + downs + [None] * (4 - len(downs))
        tables = {c: C.class_age_logcpm(C.ns(PRIMARY, c)).get(cls) for c in C.CHEMISTRIES}
        fig, axes = plt.subplots(2, 4, figsize=(12, 5), squeeze=False)
        for ax, gene in zip(axes.flat, slots):
            if gene is None:
                ax.axis("off")
                continue
            for c, marker in zip(C.CHEMISTRIES, ["o", "^"]):
                t = tables.get(c)
                if t is None or gene not in t[0].index:
                    continue
                ax.plot(t[2], t[0].loc[gene].to_numpy(), marker=marker, lw=1, ms=4, label=c)
            ax.set_title(gene, fontsize=9)
            ax.tick_params(labelsize=7)
        next(a for a, g in zip(axes.flat, slots) if g is not None).legend(fontsize=7)
        fig.suptitle(f"{PRIMARY} {cls}: top replicated rising (top row) and falling (bottom row); "
                     "y = log2 TMM-CPM, x = age (pcw)", fontsize=9)
        fig.tight_layout()
        out.figure(fig, f"trends_{PRIMARY}_{cls.replace(' ', '_')}",
                   f"Top replicated age trends in {PRIMARY} {cls}")
        plt.close(fig)


if __name__ == "__main__":
    main()
