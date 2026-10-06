#!/usr/bin/env python3
"""19 - Robustness checks: gene length, quality metrics, dissociation stress.

Question: are the main list-level results explained by three known artefacts?
  1. Gene length. Autism genes are long, and long genes behave differently in
     3' single-cell data (more intronic reads, higher in neurons). Do 06's
     class preference, 06's age coordination and 07's coherence hold when the
     random comparison genes are matched on length as well as expression?
  2. Quality metrics. Does any 08 module track mitochondrial fraction,
     unspliced fraction, UMIs or genes per cell, or doublet score across
     clusters -- within cell classes, so that marking a cell type with more
     UMIs does not count?
  3. Dissociation stress. Immediate-early and heat-shock genes rise when
     tissue is stressed during dissociation. Does a stress score change with
     age inside a cell class, and do list age trends survive controlling for it?

Method
  1. Length = genomic span (gene_id_map). Each test is run twice in the same
     code: random genes matched on expression only (as 06/07), and matched on
     expression x length tertile. Class preference: _common.set_class_preference;
     age coordination: mean rho with age (03) vs matched sets; coherence:
     across-cluster mean pairwise Spearman (07) vs matched sets. v2 x v3
     combined and tiered for both versions.
  2. Per stratum: module score per cluster = mean Z of its genes (07/08
     clusters); Spearman with the cluster's median QC values (04), raw and
     within class (class means removed from both). Flag: |rho within class|
     >= 0.5 in both chemistries, same sign.
  3. Stress score per class x age = mean Z (across the class's ages) of FOS,
     FOSB, JUN, JUNB, EGR1, IER2, IER3, ATF3, DUSP1, ZFP36, HSPA1A, HSPA1B,
     HSPA8, HSP90AA1 (van den Brink et al. 2017). Spearman with age, exact
     permutation, v2 x v3 combined. For each list: mean Spearman of its genes
     with age, plain and partial on the stress score.

Inputs:
  csv_exports/_cross_dataset/gene_id_map.csv
  csv_exports/<ds>__<chem>/09_pseudobulk/{cell_class_x_age,<finest clustering>}__pseudobulk_counts.csv
  csv_exports/<ds>__<chem>/04_clusters/cluster_profile_<clustering>.csv
  results/03_age_trends_within_cell_class/age_trends_per_stratum.csv
  results/08_coexpression_modules/modules.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "19_robustness_checks"
TITLE = "Robustness checks: gene length, quality metrics, dissociation stress"
N_RANDOM = 2000
N_RANDOM_COHERENCE = 500
MIN_GENES = 5
QC_COLS = ["frac_mito_median", "frac_unspliced_median", "total_umis_median", "n_genes_median",
           "doublet_score_median"]
QC_FLAG = 0.5
STRESS = ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "IER2", "IER3", "ATF3", "DUSP1", "ZFP36",
          "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1"]
MAX_LISTED = 8
VERSIONS = ("expression", "expression x length")


def length_tertile(ds: str) -> pd.Series:
    L = np.log10(C.gene_length(ds))
    return pd.qcut(L, 3, labels=["short", "medium", "long"]).astype(str).where(L.notna())


def list_sets(ds: str, out: C.Output) -> dict[str, list[str]]:
    return C.analysis_gene_sets(ds, groups=("ndd",), modules=False, out=out)


# ---------------------------------------------------------------------------
# 1. length-matched nulls
# ---------------------------------------------------------------------------
def length_profile(sets_by_ds: dict) -> pd.DataFrame:
    rows = []
    for ds, sets in sets_by_ds.items():
        L = C.gene_length(ds)
        expressed = C.gene_level(ds)
        ref = L.reindex(expressed.index[expressed > 1]).dropna()
        for name, genes in sets.items():
            g = L.reindex(genes).dropna()
            if len(g) < MIN_GENES:
                continue
            rows.append({"dataset": ds, "gene_set": name, "n_genes_with_length": len(g),
                         "median_length_kb": float(g.median() / 1e3),
                         "median_length_kb_expressed_genes": float(ref.median() / 1e3),
                         "median_percentile_vs_expressed": float((ref.to_numpy()[:, None] < g.to_numpy()[None, :]).mean(axis=0).mean())})
    return pd.DataFrame(rows)


def preference_both(out: C.Output, sets_by_ds: dict, rng: np.random.Generator) -> pd.DataFrame:
    parts = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv")
        lt = length_tertile(ds)
        for version, extra in zip(VERSIONS, (None, lt)):
            p = C.set_class_preference(n, ds, chem, sets_by_ds[ds], rng, n_random=N_RANDOM,
                                       label="gene_set", match_extra=extra)
            parts.append(p.assign(version=version))
    per = pd.concat(parts, ignore_index=True)
    return pd.concat([C.combine_preference(g, label="gene_set").assign(version=v)
                      for v, g in per.groupby("version")], ignore_index=True)


def age_coordination_both(out: C.Output, sets_by_ds: dict, rng: np.random.Generator) -> pd.DataFrame:
    path = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_per_stratum.csv"
    per03 = pd.read_csv(C.require(path), low_memory=False, dtype={"panels": str})
    out.used("results/03_age_trends_within_cell_class/age_trends_per_stratum.csv")
    rows = []
    for (ds, chem, cls), g in per03.groupby(["dataset", "chemistry", "cell_class"], sort=False):
        g = g.reset_index(drop=True)
        rho = g["spearman_rho"].to_numpy(float)
        dec = C.level_bins(g["mean_log2cpm"].to_numpy(float), 10)
        codes = pd.factorize(length_tertile(ds).reindex(g.gene), use_na_sentinel=True)[0]
        pos = {x: i for i, x in enumerate(g.gene)}
        for version, bins in zip(VERSIONS, (dec, dec * 5 + codes + 1)):
            for name, genes in sets_by_ds[ds].items():
                idx = np.array([pos[x] for x in genes if x in pos])
                if idx.size < MIN_GENES:
                    continue
                obs = float(rho[idx].mean())
                null = C.set_mean_rows(rho[:, None], C.matched_draws(bins, idx, N_RANDOM, rng))[:, 0]
                eff, p, _, _ = C.null_effect(obs, null)
                rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "gene_set": name,
                             "version": version, "n_ages": int(g["n_ages"].iloc[0]), "mean_rho": obs,
                             "effect_vs_null_sd": eff, "perm_p": p})
    per = pd.DataFrame(rows)
    return pd.concat([C.combine_chemistries(g, ["dataset", "cell_class", "gene_set"], weight="n_ages",
                                            labels=("rises with age", "falls with age"),
                                            carry=("mean_rho",)).assign(version=v)
                      for v, g in per.groupby("version")], ignore_index=True)


def coherence_both(out: C.Output, sets_by_ds: dict, rng: np.random.Generator, exprs: dict) -> pd.DataFrame:
    rows = []
    for ds, chem in C.STRATA:
        lc, annot = exprs[(ds, chem)]
        U = C.unit_rank_rows(lc.to_numpy(float))
        eb = C.expression_bins(lc)
        codes = pd.factorize(length_tertile(ds).reindex(lc.index), use_na_sentinel=True)[0]
        pos = {x: i for i, x in enumerate(lc.index)}
        for version, bins in zip(VERSIONS, (eb, eb * 5 + codes + 1)):
            for name, genes in sets_by_ds[ds].items():
                idx = np.array([pos[x] for x in genes if x in pos])
                if idx.size < MIN_GENES:
                    continue
                k = idx.size
                obs = float(C.coherence_from_sums(U[idx].sum(axis=0)[None, :], k)[0])
                null = C.coherence_from_sums(C.set_sums(U, C.matched_sets(bins, idx, N_RANDOM_COHERENCE, rng)), k)
                eff, p, _, _ = C.null_effect(obs, null)
                rows.append({"dataset": ds, "chemistry": chem, "gene_set": name, "version": version,
                             "coherence": obs, "effect_vs_null_sd": eff, "perm_p": p})
    per = pd.DataFrame(rows)
    return pd.concat([C.combine_chemistries(g, ["dataset", "gene_set"], labels=("coherent", "anti-coherent"),
                                            carry=("coherence",)).assign(version=v)
                      for v, g in per.groupby("version")], ignore_index=True)


def compare_versions(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    w = df.pivot_table(index=keys, columns="version", values=["stouffer_z", "tier"], aggfunc="first")
    w.columns = [f"{a}__{b}" for a, b in w.columns]
    w = w.reset_index()
    for v in VERSIONS:
        w[f"tier__{v}"] = w.get(f"tier__{v}", pd.Series("", index=w.index)).fillna("")
    w["lost_with_length"] = (w[f"tier__{VERSIONS[0]}"] != "") & (w[f"tier__{VERSIONS[1]}"] == "")
    w["gained_with_length"] = (w[f"tier__{VERSIONS[0]}"] == "") & (w[f"tier__{VERSIONS[1]}"] != "")
    return w


# ---------------------------------------------------------------------------
# 2. modules vs QC
# ---------------------------------------------------------------------------
def module_qc(out: C.Output, exprs: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = C.RESULTS / "08_coexpression_modules" / "modules.csv"
    if not path.exists():
        return pd.DataFrame(), pd.DataFrame()
    out.used("results/08_coexpression_modules/modules.csv")
    mods = pd.read_csv(path)
    mods = mods[mods.robust]
    rows = []
    for ds, chem in C.STRATA:
        lc, annot = exprs[(ds, chem)]
        Z = pd.DataFrame(C.zscore_rows(lc.to_numpy(float)), index=lc.index, columns=lc.columns)
        cls = annot.loc[lc.columns, "dominant_cell_class"].astype(str)
        prof = C.csv(C.ns(ds, chem), f"04_clusters/cluster_profile_{C.COEXPR_GROUPING[ds].replace('cluster_', '', 1)}.csv")
        prof["cluster"] = prof["cluster"].astype(str)
        prof = prof.set_index("cluster").reindex(lc.columns)
        for r in mods[mods.dataset == ds].itertuples():
            genes = [g for g in r.genes.split("|") if g in Z.index]
            if len(genes) < MIN_GENES:
                continue
            score = Z.loc[genes].mean(axis=0)
            for q in QC_COLS:
                if q not in prof or prof[q].isna().all():
                    continue
                v = pd.to_numeric(prof[q], errors="coerce")
                ok = v.notna()
                raw = float(score[ok].corr(v[ok], method="spearman"))
                sw = score[ok] - score[ok].groupby(cls[ok]).transform("mean")
                vw = v[ok] - v[ok].groupby(cls[ok]).transform("mean")
                within = float(sw.corr(vw, method="spearman"))
                rows.append({"dataset": ds, "chemistry": chem, "module": r.module, "qc_metric": q.replace("_median", ""),
                             "rho_raw": raw, "rho_within_class": within,
                             "annotation": getattr(r, "annotation", "")})
    per = pd.DataFrame(rows)
    if per.empty:
        return per, per
    w = per.pivot_table(index=["dataset", "module", "qc_metric"], columns="chemistry",
                        values=["rho_raw", "rho_within_class"]).reset_index()
    w.columns = ["_".join(c).strip("_") for c in w.columns]
    a, b = w.get("rho_within_class_v2"), w.get("rho_within_class_v3")
    w["flag_tracks_qc"] = (a.abs() >= QC_FLAG) & (b.abs() >= QC_FLAG) & (np.sign(a) == np.sign(b))
    return per, w


# ---------------------------------------------------------------------------
# 3. dissociation stress
# ---------------------------------------------------------------------------
def stress(out: C.Output, sets_by_ds: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, lrows = [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        for cls, (lc, _, ages, _) in C.class_age_logcpm(n).items():
            if len(ages) < 5:
                continue
            g = [x for x in STRESS if x in lc.index and lc.loc[x].mean() >= np.log2(5 + 1)]
            if len(g) < 3:
                continue
            Z = pd.DataFrame(C.zscore_rows(lc.to_numpy(float)), index=lc.index, columns=lc.columns)
            score = Z.loc[g].mean(axis=0).to_numpy()
            rho = float(C.spearman_rows(score[None, :], ages)[0])
            p, _ = C.spearman_perm_p(np.array([rho]), ages)
            rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "n_ages": len(ages),
                         "stress_genes": "|".join(g), "rho_vs_age": rho, "perm_p": float(p[0])})
            expressed = lc.loc[lc.mean(axis=1) >= np.log2(5 + 1)]
            plain = pd.Series(C.spearman_rows(expressed.to_numpy(float), ages), index=expressed.index)
            partial = pd.Series(C.partial_spearman(expressed.to_numpy(float), ages, score), index=expressed.index)
            for name, genes in sets_by_ds[ds].items():
                gg = [x for x in genes if x in plain.index and x not in STRESS]
                if len(gg) < MIN_GENES:
                    continue
                lrows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "gene_set": name,
                              "n_genes": len(gg), "mean_rho_age": float(plain[gg].mean()),
                              "mean_partial_rho_age_given_stress": float(partial[gg].mean())})
    per = pd.DataFrame(rows)
    comb = C.combine_chemistries(per, ["dataset", "cell_class"], effect="rho_vs_age", weight="n_ages",
                                 labels=("rises with age", "falls with age"), carry=("rho_vs_age",)) \
        if not per.empty else per
    return comb, pd.DataFrame(lrows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    out.used("_cross_dataset/gene_id_map.csv")
    sets_by_ds = {ds: list_sets(ds, out) for ds in C.DATASETS}
    exprs = {}
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        g = C.COEXPR_GROUPING[ds]
        out.used(f"{n}/09_pseudobulk/{g}__pseudobulk_counts.csv",
                 f"{n}/04_clusters/cluster_profile_{g.replace('cluster_', '', 1)}.csv")
        exprs[(ds, chem)] = C.cluster_expression(n)

    lp = length_profile(sets_by_ds)
    out.write(lp, "list_gene_length", "Per set: median gene length vs expressed genes (percentile 0.5 = typical)")
    pref = preference_both(out, sets_by_ds, rng)
    out.write(pref, "length_check_class_preference", "Class preference (as 06) with expression-only and "
              "expression x length matched nulls")
    age = age_coordination_both(out, sets_by_ds, rng)
    out.write(age, "length_check_age_coordination", "Age coordination (as 06) with both nulls")
    coh = coherence_both(out, sets_by_ds, rng, exprs)
    out.write(coh, "length_check_coherence", "Across-cluster coherence (as 07) with both nulls")
    cmp_pref = compare_versions(pref, ["dataset", "gene_set", "cell_class"])
    cmp_age = compare_versions(age, ["dataset", "cell_class", "gene_set"])
    cmp_coh = compare_versions(coh, ["dataset", "gene_set"])
    out.write(pd.concat([cmp_pref.assign(test="class preference"), cmp_age.assign(test="age coordination"),
                         cmp_coh.assign(test="coherence")], ignore_index=True), "length_check_summary",
              "Per test: tier with expression-only vs expression x length nulls; lost / gained with length")
    qc_per, qc = module_qc(out, exprs)
    if not qc.empty:
        out.write(qc_per, "module_qc_per_stratum", "Per module x QC metric x stratum: Spearman across "
                  "clusters, raw and within class")
        out.write(qc, "module_qc", f"Per module x QC metric: v2/v3 rhos; flag = |within-class rho| >= {QC_FLAG} "
                  "in both, same sign")
    st, st_lists = stress(out, sets_by_ds)
    out.write(st, "stress_vs_age", "Dissociation-stress score vs age per class; v2 x v3 combined")
    out.write(st_lists, "list_age_trends_given_stress", "Per list x class x stratum: mean rho of its genes "
              "with age, plain and partial on the stress score")

    # ---- findings -----------------------------------------------------------
    f = []
    if not lp.empty:
        lists = lp[lp.gene_set.str.startswith(("list:", "seed:"))]
        f.append("**Gene length** (median length percentile among expressed genes; 0.5 = typical): "
                 + "; ".join(f"{ds}: " + ", ".join(f"{r.gene_set} {r.median_percentile_vs_expressed:.2f}"
                                                   for r in g.sort_values("median_percentile_vs_expressed",
                                                                          ascending=False).itertuples())
                             for ds, g in lists.groupby("dataset")) + ".")
    for label, cmp_df in (("class preference", cmp_pref), ("age coordination", cmp_age), ("coherence", cmp_coh)):
        if cmp_df.empty:
            continue
        n_tier = int((cmp_df[f"tier__{VERSIONS[0]}"] != "").sum())
        lost = cmp_df[cmp_df.lost_with_length]
        gained = cmp_df[cmp_df.gained_with_length]
        desc = lambda r: f"{r.dataset} {r.gene_set}" + (f" in {r.cell_class}" if "cell_class" in cmp_df else "")
        f.append(f"**Length check, {label}**: {n_tier} tiered results with expression-matched nulls; "
                 f"{len(lost)} lost when length is matched too"
                 + (" (" + "; ".join(desc(r) for r in lost.head(MAX_LISTED).itertuples()) + ")" if len(lost) else "")
                 + f", {len(gained)} gained.")
    if not qc.empty:
        fl = qc[qc.flag_tracks_qc]
        f.append("**Modules tracking quality metrics within cell classes** (|rho| >= "
                 f"{QC_FLAG} in both chemistries): "
                 + ("; ".join(f"{r.dataset} {r.module} ~ {r.qc_metric} ({r.rho_within_class_v2:+.2f}/"
                              f"{r.rho_within_class_v3:+.2f})" for r in fl.itertuples())
                    if len(fl) else "none") + ".")
    if not st.empty:
        rep = st[st.tier != ""]
        f.append("**Dissociation-stress score vs age** (rho v2/v3): "
                 + "; ".join(f"{r.dataset} {r.cell_class} {r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f}"
                             + (f" ({r.direction}, {r.tier})" if r.tier else "") for r in st.itertuples()) + ".")
        if len(rep) and not st_lists.empty:
            m = st_lists.groupby(["dataset", "cell_class", "gene_set"])[
                ["mean_rho_age", "mean_partial_rho_age_given_stress"]].mean().reset_index()
            m = m.merge(rep[["dataset", "cell_class"]], on=["dataset", "cell_class"])
            shrink = m[(m.mean_rho_age.abs() >= 0.2)
                       & (m.mean_partial_rho_age_given_stress.abs() < 0.5 * m.mean_rho_age.abs())]
            f.append("**List age trends that shrink by half or more when stress is controlled** (classes whose "
                     "stress score changes with age; mean rho plain -> partial): "
                     + ("; ".join(f"{r.dataset} {r.gene_set} in {r.cell_class} {r.mean_rho_age:+.2f} -> "
                                  f"{r.mean_partial_rho_age_given_stress:+.2f}" for r in shrink.head(MAX_LISTED).itertuples())
                        if len(shrink) else "none") + ".")
    out.summary(
        TITLE,
        "Are the list-level results explained by gene length, by quality metrics, or by dissociation stress?",
        ["1. Length = genomic span; each test run with random genes matched on expression only and on "
         "expression x length tertile, in the same code: class preference (06), age coordination (06), "
         "across-cluster coherence (07); v2 x v3 combined and tiered for both.",
         f"2. Module score per cluster (mean Z) vs the cluster's median QC values, raw and within class; "
         f"flag |within-class rho| >= {QC_FLAG} in both chemistries.",
         f"3. Stress score = mean Z of {', '.join(STRESS)} per class x age; Spearman with age; list age "
         "trends plain and partial on the stress score."],
        f,
        ["Genomic span is a proxy for what matters in 3' data (intron content, transcript length); "
         "tertiles are coarse.",
         "Cluster QC medians summarise many cells; a module can still track QC inside clusters.",
         "The stress gene set also responds to real activity (immediate-early genes in neurons), so a "
         "stress trend is not necessarily technical."],
        ["B8: re-export with stricter per-cell QC and rerun the key analyses."])


if __name__ == "__main__":
    main()
