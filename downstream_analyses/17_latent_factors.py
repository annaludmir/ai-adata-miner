#!/usr/bin/env python3
"""17 - human_dev's 50 latent factors: what they mark, which lists they carry, how they change with age.

Question: the human_dev atlas stores 50 latent factors per cell (with gene
loadings). Which cell types and regions does each factor mark; which gene
lists concentrate among a factor's top-loading genes; and which factors
change with age -- beyond what the shifting cell-type mix explains?

Method
  A. Identity, two ways. Gene-based: the cell class whose pseudobulk
     expresses the factor's 100 top positive genes most (log2 TMM-CPM per
     class, genes Z-scored across classes) -- the factor model gives rare
     cell types extreme activity values, so this is the more readable label.
     Activity-based: the class with the highest mean activity among major
     classes (>= 2% of cells), and the top region. v2/v3 agreement of each
     factor's class profile; 08's QC correlations flag technical factors.
  B. Lists in factors: each factor's 100 top positive and 100 top negative
     genes (07_modules) as two gene sets. Overlap with each list against an
     expression-matched expectation: each list gene has the chance of being in
     the pole that genes of its total-UMI decile have (over all detected
     genes; NDD genes are highly expressed, and so are top-loading genes); the
     overlap's tail probability is exact (Poisson-binomial). BH over all tests.
  C. Age: per stratum, donor mean activity vs donor age (Spearman, exact
     permutation). Composition-adjusted: each donor's expected activity from
     its cell-class mix (sum of class share x class mean activity) is
     subtracted first, so what remains is change within cell types. v2 x v3
     combined for both versions (signed Stouffer, weights sqrt(donors)), BH,
     tiered.

Inputs (csv_exports/):
  human_dev/07_modules/{module_genes_long,module_summary}.csv
  human_dev__<chem>/08_factor_activity/factor_activity_Factors_by_{cell_class,region,donor}.csv
  human_dev__<chem>/08_factor_activity/factor_qc_correlation_Factors.csv
  human_dev__<chem>/02_composition/counts_cell_class_by_donor.csv
  human_dev__<chem>/05_confounds/crosstab_age_x_donor.csv; _cross_dataset/gene_id_map.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "17_latent_factors"
TITLE = "Latent factors (human_dev): identity, gene lists, and age trends within cell types"
DATASET = "human_dev"
MIN_OVERLAP = 3
MIN_CLASS_CELLS = 200
MAJOR_CLASS_SHARE = 0.02
MAX_LISTED = 8


def factor_no(name: str) -> int:
    return int(str(name).replace("Factors_", "").replace("Factor", ""))


def activity(n: str, by: str) -> pd.DataFrame:
    t = C.csv(n, f"08_factor_activity/factor_activity_Factors_by_{by}.csv")
    t = t[t.statistic == "mean"].copy()
    t["group_level"] = t["group_level"].astype(str)
    t = t.set_index("group_level")
    cols = [c for c in t.columns if c.startswith("Factors_")]
    out = t[cols].apply(pd.to_numeric, errors="coerce")
    out.columns = [factor_no(c) for c in cols]
    out["n_cells"] = t["n_cells"]
    return out


def pole_genes() -> dict[tuple[int, str], list[str]]:
    mg = C.csv(DATASET, "07_modules/module_genes_long.csv")
    return {(factor_no(f), "+" if d == "positive" else "-"): list(g.gene.astype(str))
            for (f, d), g in mg.groupby(["factor", "direction"])}


def poisson_binomial_sf(probs: np.ndarray, k: int) -> float:
    """P(X >= k) for X a sum of independent Bernoulli(probs)."""
    if k <= 0:
        return 1.0
    dist = np.zeros(k)          # P(X = 0 .. k-1) so far; mass beyond k-1 is dropped
    dist[0] = 1.0
    for p in probs:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= (1 - p)
    return float(max(0.0, 1.0 - dist.sum()))


def identity(out: C.Output, poles: dict) -> pd.DataFrame:
    rows = []
    prof = {}
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        out.used(f"{n}/08_factor_activity/factor_activity_Factors_by_cell_class.csv",
                 f"{n}/08_factor_activity/factor_activity_Factors_by_region.csv",
                 f"{n}/08_factor_activity/factor_qc_correlation_Factors.csv")
        out.used(f"{n}/09_pseudobulk/cell_class__pseudobulk_counts.csv")
        cl = activity(n, "cell_class")
        major = cl.index[cl.n_cells >= MAJOR_CLASS_SHARE * cl.n_cells.sum()]
        cl = cl[cl.n_cells >= MIN_CLASS_CELLS].drop(columns="n_cells")
        expr = C.tmm_log_cpm(C.group_matrix(n, "cell_class", "pseudobulk_counts", min_cells=MIN_CLASS_CELLS))
        ez = pd.DataFrame(C.zscore_rows(expr.to_numpy(float)), index=expr.index, columns=expr.columns)
        rg = activity(n, "region").drop(columns="n_cells")
        qc = C.csv(n, "08_factor_activity/factor_qc_correlation_Factors.csv")
        qc["factor"] = qc["factor"].map(factor_no)
        qc = qc.set_index("factor")
        z = (cl - cl.mean()) / cl.std().replace(0, np.nan)
        zr = (rg - rg.mean()) / rg.std().replace(0, np.nan)
        prof[chem] = z
        for k in cl.columns:
            g = [x for x in poles.get((k, "+"), []) if x in ez.index]
            gene_prof = ez.loc[g].mean(axis=0) if g else pd.Series(dtype=float)
            rows.append({"chemistry": chem, "factor": k,
                         "genes_top_cell_class": gene_prof.idxmax() if len(gene_prof) else "",
                         "genes_top_class_z": float(gene_prof.max()) if len(gene_prof) else np.nan,
                         "top_cell_class": cl.loc[cl.index.intersection(major), k].idxmax(),
                         "top_class_z": float(z[k].max()), "bottom_cell_class": z[k].idxmin(),
                         "top_region": zr[k].idxmax(), "top_region_z": float(zr[k].max()),
                         "max_abs_qc_correlation": float(qc.loc[k, "max_abs_technical_correlation"])
                         if k in qc.index else np.nan,
                         "likely_technical": bool(qc.loc[k, "likely_technical"]) if k in qc.index else False,
                         "tracks_cell_cycle": bool(qc.loc[k, "tracks_cell_cycle"]) if k in qc.index else False})
    ident = pd.DataFrame(rows)
    wide = ident.pivot(index="factor", columns="chemistry")
    res = pd.DataFrame({"factor": wide.index})
    for col in ("genes_top_cell_class", "genes_top_class_z", "top_cell_class", "top_class_z",
                "bottom_cell_class", "top_region", "max_abs_qc_correlation",
                "likely_technical", "tracks_cell_cycle"):
        for chem in C.CHEMISTRIES:
            res[f"{col}_{chem}"] = wide[(col, chem)].to_numpy()
    shared = prof["v2"].index.intersection(prof["v3"].index)
    res["class_profile_r_v2_v3"] = [float(prof["v2"].loc[shared, k].corr(prof["v3"].loc[shared, k]))
                                    for k in res.factor]
    res["same_top_class"] = res.genes_top_cell_class_v2 == res.genes_top_cell_class_v3
    return res


def lists_in_factors(out: C.Output, poles_genes: dict) -> pd.DataFrame:
    out.used(f"{DATASET}/07_modules/module_genes_long.csv", "_cross_dataset/gene_id_map.csv")
    gm = C.csv("_cross_dataset", "gene_id_map.csv", low_memory=False).query("dataset == @DATASET")
    umis = pd.to_numeric(gm["GeneTotalUMIs"], errors="coerce")
    universe = gm.loc[umis > 0].drop_duplicates("symbol")
    genes = universe["symbol"].astype(str).to_numpy()
    pos = {g: i for i, g in enumerate(genes)}
    bins = C.level_bins(np.log1p(pd.to_numeric(universe["GeneTotalUMIs"]).to_numpy(float)), 10)
    poles = {}
    for key, members in poles_genes.items():
        mask = np.zeros(len(genes), dtype=bool)
        mask[[pos[x] for x in members if x in pos]] = True
        poles[key] = mask
    sets = C.analysis_gene_sets(DATASET, groups=("ndd", "cell_cycle"), modules=False, out=out)
    rows = []
    for name, members in sets.items():
        idx = np.array(sorted({pos[x] for x in members if x in pos}))
        if idx.size < 5:
            continue
        for (k, d), mask in poles.items():
            rate = pd.Series(mask).groupby(bins).mean()          # pole share per expression decile
            probs = rate.reindex(bins[idx]).to_numpy(float)
            obs = int(mask[idx].sum())
            exp = float(probs.sum())
            rows.append({"gene_set": name, "factor": k, "pole": d, "n_set_genes": int(idx.size),
                         "pole_size": int(mask.sum()), "overlap": obs, "expected": exp,
                         "fold": obs / exp if exp > 0 else np.nan,
                         "p": poisson_binomial_sf(probs, obs),
                         "genes": "|".join(genes[idx][mask[idx]][:30])})
    res = pd.DataFrame(rows)
    res["q"] = C.bh(res.p)
    return res.sort_values("p")


def age_trends(out: C.Output) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        out.used(f"{n}/08_factor_activity/factor_activity_Factors_by_donor.csv",
                 f"{n}/02_composition/counts_cell_class_by_donor.csv", f"{n}/05_confounds/crosstab_age_x_donor.csv")
        don = C.donor_ages(n).set_index("donor")
        act = activity(n, "donor").drop(columns="n_cells")
        act = act.loc[[d for d in act.index if d in don.index]]
        cls_mean = activity(n, "cell_class").drop(columns="n_cells")
        comp = C.csv(n, "02_composition/counts_cell_class_by_donor.csv", index_col=0)
        comp.index = comp.index.astype(str)
        comp = comp.reindex(act.index)[[c for c in comp.columns if c in cls_mean.index]]
        share = comp.div(comp.sum(axis=1), axis=0)
        expected = share.to_numpy(float) @ cls_mean.loc[share.columns, act.columns].to_numpy(float)
        adjusted = act - expected
        ages = don.loc[act.index, "age_pcw"].to_numpy(float)
        for kind, M in (("raw", act), ("composition-adjusted", adjusted)):
            rho = C.spearman_rows(M.T.to_numpy(float), ages)
            p, _ = C.spearman_perm_p(rho, ages)
            for k, r, pp in zip(M.columns, rho, p):
                rows.append({"dataset": DATASET, "chemistry": chem, "version": kind, "factor": k,
                             "n_donors": len(ages), "rho_vs_age": r, "perm_p": pp})
    per = pd.DataFrame(rows)
    comb = C.combine_chemistries(per, ["dataset", "version", "factor"], effect="rho_vs_age",
                                 weight="n_donors", labels=("rises with age", "falls with age"),
                                 carry=("rho_vs_age",))
    return per, comb


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    if not (C.EXPORTS / DATASET / "07_modules" / "module_genes_long.csv").exists():
        out.summary(TITLE, "Skipped: no factor loadings exported.", [], ["**No 07_modules export.**"], [], [])
        return
    poles = pole_genes()
    ident = identity(out, poles)
    out.write(ident, "factor_identity", "Per factor: class whose pseudobulk expresses its top positive "
              "genes most; top major class by mean activity; top region; v2/v3 agreement; QC flags from 08")
    lif = lists_in_factors(out, poles)
    out.write(lif, "lists_in_factors", "Per set x factor pole: overlap with the factor's 100 top genes vs "
              "the expression-matched expectation; exact Poisson-binomial p; BH over all tests")
    per, comb = age_trends(out)
    out.write(per, "factor_age_trends_per_stratum", "Per stratum x factor: donor mean activity vs age, "
              "raw and composition-adjusted (exact permutation p)")
    out.write(comb, "factor_age_trends_combined", "v2 x v3 combined factor age trends; tier")

    lab = ident.set_index("factor")

    def describe(k: int) -> str:
        r = lab.loc[k]
        tech = " [technical?]" if r.likely_technical_v2 or r.likely_technical_v3 else ""
        cyc = " [cell cycle]" if r.tracks_cell_cycle_v2 or r.tracks_cell_cycle_v3 else ""
        return (f"F{k} ({r.genes_top_cell_class_v2}"
                f"{'' if r.same_top_class else '/' + str(r.genes_top_cell_class_v3)} genes){tech}{cyc}")

    # ---- findings -----------------------------------------------------------
    f = []
    f.append(f"**Factor identity**: {int(ident.same_top_class.sum())} of {len(ident)} factors' top genes peak in "
             f"the same cell class in both chemistries (median activity class-profile r "
             f"{ident.class_profile_r_v2_v3.median():.2f}); "
             f"{int((ident.likely_technical_v2 | ident.likely_technical_v3).sum())} flagged as likely technical "
             f"and {int((ident.tracks_cell_cycle_v2 | ident.tracks_cell_cycle_v3).sum())} as tracking the cell "
             "cycle by 08. Class where each factor's top genes peak (v2): " + ", ".join(
                 f"{c} {g.factor.size}" for c, g in ident.groupby("genes_top_cell_class_v2")) + ".")
    hits = lif[(lif.q < 0.05) & (lif.overlap >= MIN_OVERLAP)]
    if len(hits):
        parts = []
        for name, g in hits.groupby("gene_set", sort=False):
            parts.append(f"{name}: " + ", ".join(f"{describe(r.factor)}{r.pole} ({r.overlap} genes, "
                                                  f"{r.fold:.1f}x)" for r in g.head(3).itertuples()))
        f.append("**Gene sets concentrated among a factor's top genes** (q < 0.05, >= "
                 f"{MIN_OVERLAP} genes; + / - = positive / negative pole): " + " | ".join(parts[:MAX_LISTED * 2]) + ".")
    else:
        f.append("**No gene set concentrates among any factor's top genes** beyond expression-matched chance.")
    for kind in ("raw", "composition-adjusted"):
        g = comb[(comb.version == kind) & (comb.tier != "")] if not comb.empty else comb
        f.append(f"**Factors changing with age ({kind})**: "
                 + ("; ".join(f"{describe(r.factor)} {r.direction} (rho {r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f}; {r.tier})"
                              for r in g.head(MAX_LISTED).itertuples())
                    + ("" if len(g) <= MAX_LISTED else f" (+{len(g) - MAX_LISTED} more)")
                    if len(g) else "none replicated") + ".")
    if not comb.empty:
        w = comb.pivot(index="factor", columns="version", values="tier").fillna("")
        lost = w[(w["raw"] != "") & (w["composition-adjusted"] == "")].index
        f.append("**Explained by cell-type mix** (age trend in raw activity that vanishes after composition "
                 "adjustment): " + (", ".join(describe(k) for k in lost[:MAX_LISTED * 2]) if len(lost) else "none") + ".")
    out.summary(
        TITLE,
        "Which cell types and regions does each latent factor mark, which gene lists concentrate "
        "among its top genes, and which factors change with age beyond the shifting cell-type mix?",
        ["A: class whose pseudobulk (log2 TMM-CPM, Z across classes) expresses a factor's top positive "
         f"genes most; top major class (>= {MAJOR_CLASS_SHARE:.0%} of cells) by mean activity; v2/v3 "
         "agreement; 08's QC-correlation flags.",
         "B: each factor's 100 top positive and negative genes (07_modules); overlap with each set vs "
         "an expectation where each set gene has its total-UMI decile's rate of pole membership; "
         "exact Poisson-binomial p; BH over all tests.",
         "C: donor mean activity vs donor age (Spearman, exact permutation), raw and after subtracting "
         "each donor's expected activity from its class mix; v2 x v3 signed Stouffer (weights "
         "sqrt(donors)), BH, tiered."],
        f,
        ["Factor activity per donor pools all of the donor's cells; the composition adjustment uses "
         "class means from all donors, so within-class differences between donors that are not age "
         "remain in the residual.",
         "Top-gene sets are fixed at 100 genes per pole by the export; a factor whose signal is spread "
         "over many genes is under-represented.",
         "human_dev donors were dissected differently; region mix is not adjusted (only class mix)."],
        ["Export class x donor factor activity to test age within each class directly."])


if __name__ == "__main__":
    main()
