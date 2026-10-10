#!/usr/bin/env python3
"""36 - Spectra gene programmes: which hold up, what they add to the NDD lists, how they change with age.

Question: Spectra (scripts 24-25) fitted knowledge-guided gene programmes to
a cell subsample per stratum. A programme steered by a prior set can echo its
prior whatever the data say, so its labels prove nothing on their own. Which
programmes are supported by data the fit never saw; does each replicate between
the v2 and v3 fits; what did Spectra add to or drop from the NDD lists, and do
the added genes truly co-vary with the list; and how does each programme's
activity change with age inside a cell class?

Method
  A. Per stratum: factors by scope (global / cell class) and label (prior set
     recovered by its top-50 genes at overlap >= 0.2, or 'new').
  B. Replication: v2 and v3 factors paired by cosine similarity of their gene
     weights over shared genes (mutual best, >= 0.5).
  C. Held-out support: a factor's top-50 genes' mean pairwise Spearman across
     the OTHER chemistry's cluster pseudobulks (07/08 rules) -- independent
     donors -- against random gene sets matched on mean level x spread. For
     factors labelled with a prior set: the genes Spectra added (top genes not in
     the set) correlated with the set genes it kept, against random genes from
     the same bins: do the additions co-vary with the programme in other donors?
  D. Activity over age: mean cell score per class x age (script 25); Spearman
     with age per class (>= 5 ages), exact permutation; replicated factor pairs
     combined across chemistries (signed Stouffer, BH, tiered); pooled over
     regions and within each region (human_dev's sampled regions differ by age).
  E. NDD lists in data-driven programmes (top genes reproduce no prior set,
     overlap < 0.5): list genes among the top 50 vs 2,000 random sets matched on
     level x spread in the same chemistry's clusters; BH; replicated when q < 0.05
     for both factors of a v2 / v3 pair.
  E2. For list-labelled factor pairs, genes added in both fits and
     list genes kept in both (weight above the 95th percentile of non-list genes
     in the factor) and dropped in both (weight at or below their median).

Inputs: csv_exports/<ds>__<chem>/25_spectra/*; cluster pseudobulks (07/08 rules)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "36_spectra_programs"
TITLE = "Spectra gene programmes: held-out support, replication, list refinement and age trends"
N_TOP = 50
KEEP_Q = 0.95     # a list gene is kept when it outweighs 95% of non-list genes in the factor
MATCH_MIN = 0.5
N_RANDOM = 500
MIN_AGES = 5
MIN_SCORE_CELLS = 20
DATA_DRIVEN_MAX = 0.5    # top-gene overlap with every prior set below this: the factor is data-driven
MIN_SET_GENES = 10
N_RANDOM_ENRICH = 2000
MAX_LISTED = 10


def load(n: str, out: C.Output):
    base = C.EXPORTS / n / "25_spectra"
    if not (base / "factor_info.csv").exists():
        return None
    for f in ("factor_info", "factor_gene_weights", "factor_scores_by_group", "prior_sets"):
        out.used(f"{n}/25_spectra/{f}.csv")
    info = pd.read_csv(base / "factor_info.csv").set_index("factor")
    w = pd.read_csv(base / "factor_gene_weights.csv").set_index("factor")
    scores = pd.read_csv(base / "factor_scores_by_group.csv")
    pri = pd.read_csv(base / "prior_sets.csv")
    return info, w, scores, pri


def pair_factors(w2: pd.DataFrame, w3: pd.DataFrame) -> pd.DataFrame:
    genes = w2.columns.intersection(w3.columns)
    a, b = w2[genes].to_numpy(float), w3[genes].to_numpy(float)
    a = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
    b = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-12)
    S = C.dot(a, b.T)
    best3, best2 = S.argmax(axis=1), S.argmax(axis=0)
    rows = []
    for i, j in enumerate(best3):
        if best2[j] == i and S[i, j] >= MATCH_MIN:
            rows.append({"factor_v2": w2.index[i], "factor_v3": w3.index[j], "cosine": float(S[i, j])})
    return pd.DataFrame(rows)


def held_out(info: pd.DataFrame, w: pd.DataFrame, prior_genes: dict, lc: pd.DataFrame,
             rng: np.random.Generator) -> pd.DataFrame:
    U = C.unit_rank_rows(lc.to_numpy(float))
    bins = C.expression_bins(lc)
    pos = {g: i for i, g in enumerate(lc.index)}
    rows = []
    for f, r in info.iterrows():
        top = list(w.loc[f].sort_values(ascending=False).index[:N_TOP])
        idx = np.array([pos[g] for g in top if g in pos])
        row = {"factor": f, "label": r.label, "n_top_in_clusters": int(idx.size)}
        if idx.size >= 5:
            k = idx.size
            coh = float(C.coherence_from_sums(U[idx].sum(axis=0)[None, :], k)[0])
            null = C.coherence_from_sums(C.set_sums(U, C.matched_sets(bins, idx, N_RANDOM, rng)), k)
            eff, p, mu, _ = C.null_effect(coh, null)
            row.update({"coherence": coh, "null_mean": mu, "effect_vs_null_sd": eff, "perm_p": p})
        prior = prior_genes.get(r.label)
        if prior:
            kept = [g for g in top if g in prior and g in pos]
            added = [g for g in top if g not in prior and g in pos]
            if len(kept) >= 3 and len(added) >= 3:
                ki, ai = np.array([pos[g] for g in kept]), np.array([pos[g] for g in added])
                obs = float((U[ai] @ U[ki].sum(axis=0)).mean() / len(ki))
                rand = C.matched_sets(bins, ai, N_RANDOM, rng)
                null = (C.set_sums(U, rand) @ U[ki].sum(axis=0)) / (len(ki) * ai.size)
                eff, p, mu, _ = C.null_effect(obs, null)
                row.update({"n_kept": len(kept), "n_added": len(added), "added_vs_kept_r": obs,
                            "added_null_mean": mu, "added_effect_vs_null_sd": eff, "added_perm_p": p,
                            "added_genes": "|".join(added[:25])})
        rows.append(row)
    return pd.DataFrame(rows)


def age_trends(scores: pd.DataFrame, ds: str) -> pd.DataFrame:
    """Spearman of mean cell score with age per factor x class, pooled over regions and within each region."""
    parts = []
    s = scores[scores.grouping == "cell_class_x_age"].copy()
    sp = s.group.str.split(r" \| ", regex=True)
    s["cell_class"], s["region"], s["age"] = sp.str[0], "all regions", sp.str[1]
    parts.append(s)
    r = scores[scores.grouping == "cell_class_x_region_x_age"].copy()
    if len(r):
        sp = r.group.str.split(r" \| ", regex=True)
        r["cell_class"], r["region"], r["age"] = sp.str[0], sp.str[1], sp.str[2]
        parts.append(r)
    s = pd.concat(parts, ignore_index=True)
    s["age"] = pd.to_numeric(s.age, errors="coerce")
    s = s[(s.n_cells >= MIN_SCORE_CELLS) & s.age.notna()]
    s = s[~s.age.map(lambda a: C.excluded(ds, age=a))]
    rows = []
    for (cls, reg), g in s.groupby(["cell_class", "region"]):
        M = g.pivot_table(index="factor", columns="age", values="mean_score").dropna(axis=1)
        if M.shape[1] < MIN_AGES:
            continue
        x = M.columns.to_numpy(float)
        rho = C.spearman_rows(M.to_numpy(float), x)
        p, _ = C.spearman_perm_p(rho, x)
        rows.append(pd.DataFrame({"factor": M.index, "cell_class": cls, "region": reg, "n_ages": M.shape[1],
                                  "rho_vs_age": rho, "perm_p": p}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def ndd_enrichment(info: pd.DataFrame, w: pd.DataFrame, sets: dict, lc: pd.DataFrame,
                   rng: np.random.Generator) -> pd.DataFrame:
    """NDD sets among the top genes of data-driven factors, against expression-matched random sets.

    Data-driven = the factor's top genes do not reproduce a prior set (overlap
    coefficient < DATA_DRIVEN_MAX). NDD genes are long and highly expressed in
    neurons, so the null draws random sets matched on level x spread across the
    same chemistry's clusters; genes and background = model genes in those clusters.
    """
    genes = [g for g in w.columns if g in lc.index]
    pos = {g: i for i, g in enumerate(genes)}
    bins = C.expression_bins(lc.loc[genes])
    dd = info[info.overlap_coefficient < DATA_DRIVEN_MAX]
    rows = []
    for name, members in sets.items():
        idx = np.array(sorted(pos[g] for g in members if g in pos))
        if idx.size < MIN_SET_GENES:
            continue
        rand = C.matched_sets(bins, idx, N_RANDOM_ENRICH, rng)
        for f in dd.index:
            top = w.loc[f, genes].nlargest(N_TOP).index
            mask = np.zeros(len(genes), bool)
            mask[[pos[g] for g in top]] = True
            k = int(mask[idx].sum())
            null = mask[rand].sum(axis=1)
            rows.append({"factor": f, "scope": info.loc[f, "scope"], "label": info.loc[f, "label"],
                         "gene_set": name, "set_genes": idx.size, "in_top": k, "null_mean": float(null.mean()),
                         "fold": k / max(float(null.mean()), 1e-9),
                         "perm_p": float((np.sum(null >= k) + 1) / (N_RANDOM_ENRICH + 1)),
                         "genes": "|".join(g for g in top if g in members)})
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    fits = {(ds, chem): load(C.ns(ds, chem), out) for ds, chem in C.STRATA}
    if all(v is None for v in fits.values()):
        out.summary(TITLE, "Skipped: no Spectra fits yet.",
                    ["Needs csv_exports/<ns>/25_spectra/ from scripts 24-25 (running_scripts/slurm_spectra.sh)."],
                    ["**Not run**: fit Spectra first (slurm_spectra.sh per dataset x chemistry), then rerun step 3."],
                    [], [])
        return
    overview, pairs_all, support, ages, refine, enrich = [], [], [], [], [], []
    tops = {}
    for ds in C.DATASETS:
        f2, f3 = fits.get((ds, "v2")), fits.get((ds, "v3"))
        for chem, fit in (("v2", f2), ("v3", f3)):
            if fit is None:
                continue
            info, w, scores, pri = fit
            overview.append(info.reset_index().assign(dataset=ds, chemistry=chem))
            for fac in w.index:
                tops[(ds, chem, fac)] = set(w.loc[fac].sort_values(ascending=False).index[:N_TOP])
            other = "v3" if chem == "v2" else "v2"
            lc, _ = C.cluster_expression(C.ns(ds, other))
            out.used(f"{C.ns(ds, other)}/09_pseudobulk/{C.COEXPR_GROUPING[ds]}__pseudobulk_counts.csv")
            # prior set genes from the same sources the fit used (labels: list:, seed:, marker:)
            # (lists uncollapsed: the fit's prior keeps every gene of a locus)
            prior_genes = {f"list:{k}": set(v) for k, v in C.mapped_lists(ds, collapse=False).items()}
            pan = C.panels(C.ns(ds, "v2"))
            out.used(f"{C.ns(ds, 'v2')}/11_panels/panel_coverage.csv")
            for (grp, pname), g in pan[pan.panel_group.isin(["ndd", "cell_cycle", "marker"])].groupby(
                    ["panel_group", "panel"]):
                prior_genes[f"{'marker' if grp == 'marker' else 'seed'}:{pname}"] = set(g.gene)
            h = held_out(info, w, prior_genes, lc, rng)
            support.append(h.assign(dataset=ds, chemistry=chem, validated_on=other))
            a = age_trends(scores, ds)
            if len(a):
                ages.append(a.assign(dataset=ds, chemistry=chem))
            own, _ = C.cluster_expression(C.ns(ds, chem))
            out.used(f"{C.ns(ds, chem)}/09_pseudobulk/{C.COEXPR_GROUPING[ds]}__pseudobulk_counts.csv")
            ndd_sets = {k: v for k, v in prior_genes.items() if k.startswith("list:")}
            for pname, g in pan[pan.panel_group == "ndd"].groupby("panel"):
                ndd_sets[f"seed:{pname}"] = set(g.gene)
            e = ndd_enrichment(info, w, ndd_sets, own, rng)
            if len(e):
                enrich.append(e.assign(dataset=ds, chemistry=chem))
        if f2 is not None and f3 is not None:
            p = pair_factors(f2[1], f3[1])
            p["label_v2"] = p.factor_v2.map(f2[0].label)
            p["label_v3"] = p.factor_v3.map(f3[0].label)
            p["scope_v2"] = p.factor_v2.map(f2[0].scope)
            p["scope_v3"] = p.factor_v3.map(f3[0].scope)
            pairs_all.append(p.assign(dataset=ds))
            for r in p.itertuples():
                if str(r.label_v2).startswith("list:") and r.label_v2 == r.label_v3:
                    t2 = set(f2[1].loc[r.factor_v2].sort_values(ascending=False).index[:N_TOP])
                    t3 = set(f3[1].loc[r.factor_v3].sort_values(ascending=False).index[:N_TOP])
                    members = set(C.mapped_lists(ds, collapse=False).get(r.label_v2.split(":", 1)[1], []))
                    in_model = members & set(f2[1].columns) & set(f3[1].columns)
                    if not in_model:          # list no longer among the gene lists
                        continue
                    kept, dropped = [], []
                    for fit, fac in ((f2, r.factor_v2), (f3, r.factor_v3)):
                        wt = fit[1].loc[fac]
                        bg = wt[~wt.index.isin(list(members))]
                        kept.append({g for g in in_model if wt[g] > bg.quantile(KEEP_Q)})
                        dropped.append({g for g in in_model if wt[g] <= bg.median()})
                    k_both, d_both = kept[0] & kept[1], dropped[0] & dropped[1]
                    refine.append({"dataset": ds, "gene_list": r.label_v2, "factor_v2": r.factor_v2,
                                   "factor_v3": r.factor_v3, "cosine": r.cosine,
                                   "list_genes_in_model": len(in_model),
                                   "added_in_both": "|".join(sorted((t2 & t3) - members)),
                                   "kept_in_both": "|".join(sorted(k_both, key=lambda g: -f2[1].loc[r.factor_v2, g])),
                                   "dropped_in_both": "|".join(sorted(d_both)),
                                   "n_added_in_both": len((t2 & t3) - members),
                                   "n_kept_in_both": len(k_both), "n_dropped_in_both": len(d_both)})
    overview = pd.concat(overview, ignore_index=True)
    pairs = pd.concat(pairs_all, ignore_index=True) if pairs_all else pd.DataFrame()
    support = pd.concat(support, ignore_index=True)
    ages = pd.concat(ages, ignore_index=True) if ages else pd.DataFrame()
    refine = pd.DataFrame(refine)
    enrich = pd.concat(enrich, ignore_index=True) if enrich else pd.DataFrame()
    erep = pd.DataFrame()
    if len(enrich):
        enrich["q"] = np.nan
        for _, ix in enrich.groupby("dataset").groups.items():
            enrich.loc[ix, "q"] = C.bh(enrich.loc[ix, "perm_p"])
        rows = []
        for p in pairs.itertuples():
            e2 = enrich[(enrich.dataset == p.dataset) & (enrich.chemistry == "v2") & (enrich.factor == p.factor_v2)]
            e3 = enrich[(enrich.dataset == p.dataset) & (enrich.chemistry == "v3") & (enrich.factor == p.factor_v3)]
            m = e2.merge(e3, on="gene_set", suffixes=("_v2", "_v3"))
            for r in m.itertuples():
                rows.append({"dataset": p.dataset, "pair": f"{p.factor_v2}~{p.factor_v3}", "scope": p.scope_v2,
                             "cosine": p.cosine, "gene_set": r.gene_set, "fold_v2": r.fold_v2, "fold_v3": r.fold_v3,
                             "q_v2": r.q_v2, "q_v3": r.q_v3,
                             "replicated": bool(r.q_v2 < 0.05 and r.q_v3 < 0.05),
                             "genes_in_both": "|".join(sorted(set(r.genes_v2.split("|")) & set(r.genes_v3.split("|"))
                                                              - {""}))})
        erep = pd.DataFrame(rows)
    acomb = pd.DataFrame()
    if len(pairs) and len(ages):
        rows = []
        for p in pairs.itertuples():
            for chem, f in (("v2", p.factor_v2), ("v3", p.factor_v3)):
                a = ages[(ages.dataset == p.dataset) & (ages.chemistry == chem) & (ages.factor == f)]
                for r in a.itertuples():
                    rows.append({"dataset": p.dataset, "pair": f"{p.factor_v2}~{p.factor_v3}",
                                 "label": p.label_v2 if p.label_v2 == p.label_v3 else f"{p.label_v2} / {p.label_v3}",
                                 "scope": p.scope_v2, "cell_class": r.cell_class, "region": r.region,
                                 "chemistry": chem, "rho_vs_age": r.rho_vs_age,
                                 "perm_p": r.perm_p, "n_ages": r.n_ages})
        if rows:
            acomb = C.combine_chemistries(pd.DataFrame(rows), ["dataset", "pair", "label", "scope", "cell_class",
                                                               "region"],
                                          effect="rho_vs_age", weight="n_ages",
                                          labels=("rises with age", "falls with age"), carry=("rho_vs_age",))
    support["q"] = C.bh(support.perm_p) if "perm_p" in support else np.nan
    if "added_perm_p" in support:
        support["added_q"] = C.bh(support.added_perm_p)
    mod_path = C.RESULTS / "08_coexpression_modules" / "modules.csv"
    if mod_path.exists():
        out.used("results/08_coexpression_modules/modules.csv")
        mods = pd.read_csv(mod_path)
        mods = mods[mods.robust]
        best = []
        for r in overview.itertuples():
            top = tops[(r.dataset, r.chemistry, r.factor)]
            cand = [(len(top & set(m.genes.split("|"))) / max(1, min(len(top), len(m.genes.split("|")))), m.module)
                    for m in mods[mods.dataset == r.dataset].itertuples()]
            best.append(max(cand) if cand else (np.nan, ""))
        overview["best_08_module"] = [b[1] for b in best]
        overview["best_08_module_overlap"] = [b[0] for b in best]
    out.write(overview, "factors", "Per stratum x factor: scope, label, overlap with its prior set, top genes")
    out.write(pairs, "factor_pairs_v2_v3", f"Mutual-best v2 / v3 factor pairs by cosine of gene weights (>= {MATCH_MIN})")
    out.write(support, "held_out_support", "Per factor: top-gene coherence across the other chemistry's clusters vs "
              "matched genes; added-vs-kept genes correlation for prior-labelled factors")
    out.write(ages, "factor_age_trends_per_stratum", "Per factor x class: Spearman of mean cell score with age")
    out.write(acomb, "factor_age_trends_combined", "Per replicated factor pair x class: v2 x v3 combined; tier")
    out.write(enrich, "ndd_enrichment_data_driven", "Per data-driven factor x NDD set: set genes among the top "
              f"{N_TOP} vs expression-matched random sets")
    out.write(erep, "ndd_enrichment_replicated", "Per v2 / v3 pair of data-driven factors x NDD set: both chemistries")
    out.write(refine, "list_refinement", "Per list-labelled factor pair: genes added in both fits, list genes kept / dropped in both")

    f = []
    trained = {}
    for ds, chem in C.STRATA:
        fs = C.EXPORTS / C.ns(ds, chem) / "25_spectra" / "fit_summary.csv"
        if fs.exists():
            out.used(f"{C.ns(ds, chem)}/25_spectra/fit_summary.csv")
            r = pd.read_csv(fs).iloc[0]
            if pd.notna(r.epochs_run):
                trained[(ds, chem)] = (f", {int(r.epochs_run)} of {int(r.epochs_max)} epochs"
                                       + (f" ({r.stop_reason})" if "stop_reason" in r else ""))
    ov = overview.groupby(["dataset", "chemistry"]).apply(
        lambda g: f"{len(g)} factors, {int((g.label != 'new').sum())} recovering a prior set, "
                  f"{int((g.label == 'new').sum())} new", include_groups=False)
    f.append("**Fits**: " + "; ".join(f"{ds} {ch}: {v}{trained.get((ds, ch), '')}" for (ds, ch), v in ov.items()) + ".")
    if len(pairs):
        for ds, g in pairs.groupby("dataset"):
            same = g[g.label_v2 == g.label_v3]
            new = g[(g.label_v2 == "new") & (g.label_v3 == "new")]
            f.append(f"**{ds}: v2 / v3 replication** -- {len(g)} factor pairs (cosine >= {MATCH_MIN}); {len(same)} with the "
                     f"same label; {len(new)} new (data-driven) programmes found in both fits"
                     + (": " + "; ".join(f"{r.factor_v2}~{r.factor_v3} ({r.scope_v2}, cos {r.cosine:.2f})"
                                         for r in new.head(6).itertuples()) if len(new) else "") + ".")
    if "effect_vs_null_sd" in support:
        s = support.dropna(subset=["effect_vs_null_sd"])
        for ds, g in s.groupby("dataset"):
            ok = g[(g.q < 0.05) & (g.effect_vs_null_sd > 0)]
            f.append(f"**{ds}: programmes whose top genes co-express in the other chemistry's donors** "
                     f"(q < 0.05): {len(ok)} of {len(g)} ({int((ok.label != 'new').sum())} of "
                     f"{int((g.label != 'new').sum())} prior-labelled, {int((ok.label == 'new').sum())} of "
                     f"{int((g.label == 'new').sum())} new). Supported, strongest first: " + ("; ".join(
                         f"{r.chemistry} {r.factor} {r.label} ({r.effect_vs_null_sd:+.1f} SD)"
                         for r in ok.sort_values("effect_vs_null_sd", ascending=False).head(MAX_LISTED * 2).itertuples())
                         or "none") + ".")
    if "added_effect_vs_null_sd" in support:
        a = support.dropna(subset=["added_effect_vs_null_sd"])
        a = a[a.label.astype(str).str.startswith(("list:", "seed:"))]
        if len(a):
            good = a[(a.added_q < 0.05) & (a.added_effect_vs_null_sd > 0)]
            f.append("**Genes Spectra added to lists and panels that co-vary with the kept genes in held-out donors**: "
                     + ("; ".join(f"{r.dataset} {r.chemistry} {r.label}: {int(r.n_added)} added, r {r.added_vs_kept_r:.2f} vs "
                                  f"{r.added_null_mean:.2f} ({', '.join(str(r.added_genes).split('|')[:8])})"
                                  for r in good.head(MAX_LISTED).itertuples()) if len(good) else "none") + ".")
    if len(refine):
        kept = refine.n_kept_in_both.sum() / max(refine.list_genes_in_model.sum(), 1)
        lf = overview[overview.label.astype(str).str.startswith("list:")]
        f.append(f"**List programmes reproduce their lists**: {int((lf.overlap_coefficient >= 0.99).sum())} of {len(lf)} "
                 f"list-labelled factors have top-{N_TOP} genes made only of list genes, and across replicated list "
                 f"programmes {kept:.0%} of list genes outweigh 95% of non-list genes in both fits. At lam = 0.01 the "
                 "prior holds the whole list together, so Spectra does not prune the lists; whether a list coheres in "
                 "the data is what the held-out test above (and 07) measures (list_refinement.csv has the details).")
    if len(erep):
        rp = erep[erep.replicated].copy()
        rp["fmin"] = rp[["fold_v2", "fold_v3"]].min(axis=1)
        rp = rp.sort_values("fmin", ascending=False)
        n_dd = int((overview.overlap_coefficient < DATA_DRIVEN_MAX).sum())
        v2 = overview[overview.chemistry == "v2"]
        lead = {(r.dataset, r.factor): ", ".join(str(r.top_genes).split("|")[:4]) for r in v2.itertuples()}
        f.append(f"**NDD genes in data-driven programmes** ({n_dd} factors whose top genes reproduce no prior set; "
                 f"NDD set among the top {N_TOP} vs expression-matched random sets; q < 0.05 in both chemistries): "
                 + ("; ".join(f"{r.dataset} {r.pair} ({r.scope}; top {lead[(r.dataset, r.pair.split('~')[0])]}) "
                              f"{r.gene_set} x{r.fold_v2:.1f}/x{r.fold_v3:.1f} "
                              f"({', '.join(str(r.genes_in_both).split('|')[:8])})"
                              for r in rp.head(MAX_LISTED * 2).itertuples())
                    if len(rp) else "none") + ".")
    if "best_08_module_overlap" in overview:
        nm = overview[(overview.label == "new") & (overview.best_08_module_overlap < 0.2)]
        f.append(f"**New programmes vs the 08 co-expression modules**: {int((overview.label == 'new').sum())} new "
                 f"factors, {len(nm)} overlapping no robust 08 module (top-gene overlap < 0.2)"
                 + (": " + "; ".join(f"{r.dataset} {r.chemistry} {r.factor} "
                                     f"({', '.join(str(r.top_genes).split('|')[:6])})"
                                     for r in nm.head(MAX_LISTED).itertuples()) if len(nm) else "") + ".")
    if len(acomb):
        rep = acomb[acomb.tier != ""]
        reg = rep[rep.region != "all regions"]
        f.append("**Programme activity changing with age within one region** (removes the region mix that differs by "
                 "age; rho v2/v3): " + ("; ".join(
            f"{r.dataset} {r.label} [{r.pair}, {r.scope}] in {r.region} {r.cell_class} {r.direction} "
            f"({r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f}; {r.tier})"
            for r in reg.head(MAX_LISTED * 2).itertuples()) if len(reg) else "none replicated") + ".")
        pooled = rep[rep.region == "all regions"]
        if len(pooled):
            inreg = set(zip(reg.dataset, reg.pair, reg.cell_class, reg.direction))
            pooled = pooled.assign(confirmed=[(r.dataset, r.pair, r.cell_class, r.direction) in inreg
                                              for r in pooled.itertuples()])
            f.append(f"**Pooled over regions**: {len(pooled)} replicated or supported trends, "
                     f"{int(pooled.confirmed.sum())} of them also within a region. Pooled only (may follow which regions "
                     "were sampled at each age): " + ("; ".join(
                         f"{r.dataset} {r.label} [{r.pair}] in {r.cell_class} {r.direction}"
                         for r in pooled[~pooled.confirmed].head(MAX_LISTED).itertuples()) or "none") + ".")
    out.summary(
        TITLE,
        "Which Spectra programmes are supported by donors the fit never saw, which replicate between the v2 and v3 "
        "fits, what did Spectra add to or drop from the NDD lists, and how does programme activity change with age?",
        ["Fits: scripts 24-25 (stratified ~25k-cell subsample per stratum; prior = gene lists, seed NDD and cell-cycle "
         "panels, core GO processes, class markers; extra free factors).",
         f"Replication: mutual-best cosine of gene weights >= {MATCH_MIN}. Held-out support: mean pairwise Spearman of a "
         f"factor's top {N_TOP} genes across the other chemistry's cluster pseudobulks vs {N_RANDOM} random sets matched "
         "on level x spread; added genes' correlation with kept prior genes vs matched random genes; BH.",
         f"NDD enrichment: in factors whose top genes reproduce no prior set (overlap < {DATA_DRIVEN_MAX}), NDD set "
         f"genes among the top {N_TOP} vs {N_RANDOM_ENRICH:,} random sets matched on level x spread; BH per dataset; "
         "replicated when q < 0.05 for both factors of a v2 / v3 pair.",
         f"Age: Spearman of mean cell score per class x age (>= {MIN_AGES} ages, >= {MIN_SCORE_CELLS} cells), pooled "
         "and within each region; replicated pairs combined (signed Stouffer, BH, tiered). Lists: added = top-50 in both fits and not in the list; kept = weight above "
         f"the {KEEP_Q:.0%} quantile of non-list genes in both; dropped = at or below their median in both."],
        f,
        ["A prior-steered factor can echo its prior whatever the data; only the held-out tests and the v2 / v3 "
         "replication speak to support.",
         "Fits use a subsample (~25k cells per stratum), so rare cell types contribute few cells, and per-region "
         "age points can rest on a few dozen cells.",
         "The gpu backend (Spectra_gpu, minibatched) is marked by its authors as in development; the cpu backend is "
         "the published model."],
        ["Compare with scHPF or cNMF (data-driven programmes) and with STRING-expanded priors."])


if __name__ == "__main__":
    main()
