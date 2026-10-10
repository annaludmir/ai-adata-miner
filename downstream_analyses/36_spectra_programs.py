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
  E. NDD lists vs programmes they did not shape: a second fit (SPECTRA_PRIOR=
     core: GO processes, cell-cycle panels, class markers, 15 free global
     factors; no NDD set) -- in the full fit the big lists dominate the global
     gene-gene graph that every global factor, free ones too, is fitted to. Per
     core factor: NDD set genes among its top 50 vs the expectation from each
     gene's level x spread bin in the same chemistry's clusters (exact
     Poisson-binomial p); a steered factor's own prior-set genes are left out;
     BH per dataset; replicated when q < 0.05 for both factors of a v2 / v3 pair;
     the pair's held-out support and age trend are reported with it. The same
     for a third fit (SPECTRA_PRIOR=string: core prior + STRING network modules;
     a module-steered factor's module genes are left out), plus each STRING
     module's own NDD content (vs the expression-matched expectation) next to
     whether its programme replicates v2 / v3 and holds in held-out donors.
  E2. For list-labelled factor pairs, genes added in both fits and
     list genes kept in both (weight above the 95th percentile of non-list genes
     in the factor) and dropped in both (weight at or below their median).

Inputs: csv_exports/<ds>__<chem>/25_spectra/*, 25_spectra_core/*, 25_spectra_string/*; cluster pseudobulks (07/08 rules);
GO BP gene sets (annotations)
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
N_FREE_GLOBAL = 5        # free factors per scope, as script 25 sets them (N_NEW_GLOBAL; +1 per class)
N_FREE_CLASS = 1
MIN_SET_GENES = 10
# the fits whose prior holds no NDD set: (export subfolder, output-file prefix, name in the summary)
VARIANTS = {"core": ("25_spectra_core", "core", "Core-prior fit"),
            "string": ("25_spectra_string", "string", "STRING-prior fit")}
MAX_LISTED = 10


def load(n: str, out: C.Output, sub: str = "25_spectra"):
    base = C.EXPORTS / n / sub
    if not (base / "factor_info.csv").exists():
        return None
    for f in ("factor_info", "factor_gene_weights", "factor_scores_by_group", "prior_sets"):
        out.used(f"{n}/{sub}/{f}.csv")
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


def free_factors(info: pd.DataFrame) -> pd.Series:
    """True for the free (prior-less) factors of a fit.

    Script 25 gives each scope its prior sets plus free factors -- N_FREE_GLOBAL
    after the global sets, one after each class's marker sets -- and Spectra keeps
    them in that order, so the free factors are the last of each scope. Steered
    factors that drifted from their set are not free: their NDD content may come
    from the set that seeded them.
    """
    if "free" in info:                         # written by script 25 since the core prior
        return info["free"].astype(bool)
    free = pd.Series(False, index=info.index)
    for scope, g in info.groupby("scope", sort=False):
        n = N_FREE_GLOBAL if scope == "global" else N_FREE_CLASS
        free[g.index[-n:]] = True
    return free


def ndd_enrichment(info: pd.DataFrame, w: pd.DataFrame, sets: dict, lc: pd.DataFrame,
                   steering: dict | None = None) -> pd.DataFrame:
    """NDD sets among the top genes of every factor, against expression-matched expectation.

    Meant for the core-prior fit, whose prior holds no NDD set, so no factor was
    steered towards one (in the full fit the lists dominate the global graph that
    every global factor, free ones included, is fitted to). A steered factor
    shares genes with its own prior set (the synapse GO set holds synaptic NDD
    genes, say), so for a factor that is not free, the genes of the prior set it
    matches (steering[label]) are left out of the NDD set and the background.

    NDD genes are long and highly expressed in neurons, so each set gene's chance
    of being a top gene is the top-gene share of its level x spread bin across the
    same chemistry's clusters; the overlap's tail probability is exact
    (Poisson-binomial). Genes and background = model genes in those clusters.
    """
    genes = [g for g in w.columns if g in lc.index]
    pos = {g: i for i, g in enumerate(genes)}
    bins = C.expression_bins(lc.loc[genes])
    rows = []
    free = free_factors(info)
    steering = steering or {}
    for f in info.index:
        top = w.loc[f, genes].nlargest(N_TOP).index
        mask = np.zeros(len(genes), bool)
        mask[[pos[g] for g in top]] = True
        drop = set() if free[f] else {g.upper() for g in steering.get(info.loc[f, "label"], ())}
        keep = np.array([g.upper() not in drop for g in genes])
        rate = pd.Series(mask[keep]).groupby(bins[keep]).mean()
        for name, members in sets.items():
            idx = np.array(sorted(pos[g] for g in members if g in pos and g.upper() not in drop))
            if idx.size < MIN_SET_GENES:
                continue
            probs = rate.reindex(bins[idx]).to_numpy(float)
            k, exp = int(mask[idx].sum()), float(probs.sum())
            rows.append({"factor": f, "scope": info.loc[f, "scope"], "label": info.loc[f, "label"],
                         "free": bool(free[f]), "gene_set": name, "prior_genes_left_out": len(drop), "set_genes": idx.size, "in_top": k, "expected": exp,
                         "fold": k / exp if exp > 0 else np.nan, "p": C.poisson_binomial_sf(probs, k),
                         "genes": "|".join(g for g in top if g in members)})
    return pd.DataFrame(rows)


def labelled_pairs(f2, f3, ds: str) -> pd.DataFrame:
    p = pair_factors(f2[1], f3[1])
    for chem, fit in (("v2", f2), ("v3", f3)):
        p[f"label_{chem}"] = p[f"factor_{chem}"].map(fit[0].label)
        p[f"scope_{chem}"] = p[f"factor_{chem}"].map(fit[0].scope)
    return p.assign(dataset=ds)


def ndd_sets_for(ds: str, out: C.Output) -> dict:
    """User lists (uncollapsed) and seed NDD panels."""
    sets = {f"list:{k}": set(v) for k, v in C.mapped_lists(ds, collapse=False).items()}
    pan = C.panels(C.ns(ds, "v2"))
    out.used(f"{C.ns(ds, 'v2')}/11_panels/panel_coverage.csv")
    for pname, g in pan[pan.panel_group == "ndd"].groupby("panel"):
        sets[f"seed:{pname}"] = set(g.gene)
    return sets


def steering_sets(ds: str, out: C.Output) -> dict:
    """The core prior's sets by label (go:, seed: cell cycle, marker:), as script 25 names them."""
    sets = {}
    gmt = C.annotation("gmt/GO_Biological_Process_2023.gmt")
    if gmt is not None:
        out.used("annotations/gmt/GO_Biological_Process_2023.gmt")
        for term, genes in C.read_gmt(gmt).items():
            sets["go:" + term.split(" (GO:")[0]] = set(genes)
    pan = C.panels(C.ns(ds, "v2"))
    for (grp, pname), g in pan[pan.panel_group.isin(["cell_cycle", "marker"])].groupby(["panel_group", "panel"]):
        sets[f"{'marker' if grp == 'marker' else 'seed'}:{pname}"] = set(g.gene)
    return sets


def combine_ages(pairs: pd.DataFrame, ages: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for p in pairs.itertuples():
        for chem, f in (("v2", p.factor_v2), ("v3", p.factor_v3)):
            a = ages[(ages.dataset == p.dataset) & (ages.chemistry == chem) & (ages.factor == f)]
            for r in a.itertuples():
                rows.append({"dataset": p.dataset, "pair": f"{p.factor_v2}~{p.factor_v3}",
                             "label": p.label_v2 if p.label_v2 == p.label_v3 else f"{p.label_v2} / {p.label_v3}",
                             "scope": p.scope_v2, "cell_class": r.cell_class, "region": r.region,
                             "chemistry": chem, "rho_vs_age": r.rho_vs_age, "perm_p": r.perm_p, "n_ages": r.n_ages})
    if not rows:
        return pd.DataFrame()
    return C.combine_chemistries(pd.DataFrame(rows), ["dataset", "pair", "label", "scope", "cell_class", "region"],
                                 effect="rho_vs_age", weight="n_ages", labels=("rises with age", "falls with age"),
                                 carry=("rho_vs_age",))


def string_module_sets(ds: str, sub: str, out: C.Output) -> tuple[dict, dict]:
    """STRING modules of a dataset's fits: {name: genes (union over strata)}, {name: nearest GO term}."""
    genes, go = {}, {}
    for chem in C.CHEMISTRIES:
        p = C.EXPORTS / C.ns(ds, chem) / sub / "string_modules.csv"
        if p.exists():
            out.used(f"{C.ns(ds, chem)}/{sub}/string_modules.csv")
            for r in pd.read_csv(p).itertuples():
                genes.setdefault(r.module, set()).update(str(r.genes).split("|"))
                go[r.module] = r.best_go_term
    return genes, go


def module_ndd_content(modules: dict, sets: dict, lc: pd.DataFrame) -> pd.DataFrame:
    """NDD set genes in each STRING module vs the expression-matched expectation (exact p)."""
    genes = list(lc.index)
    pos = {g: i for i, g in enumerate(genes)}
    bins = C.expression_bins(lc)
    rows = []
    for m, mg in modules.items():
        mask = np.zeros(len(genes), bool)
        mask[[pos[g] for g in mg if g in pos]] = True
        if mask.sum() < MIN_SET_GENES:
            continue
        rate = pd.Series(mask).groupby(bins).mean()
        for name, members in sets.items():
            idx = np.array(sorted(pos[g] for g in members if g in pos))
            if idx.size < MIN_SET_GENES:
                continue
            probs = rate.reindex(bins[idx]).to_numpy(float)
            k, exp = int(mask[idx].sum()), float(probs.sum())
            rows.append({"module": m, "module_genes": int(mask.sum()), "gene_set": name, "in_module": k,
                         "expected": exp, "fold": k / exp if exp > 0 else np.nan, "p": C.poisson_binomial_sf(probs, k),
                         "genes": "|".join(genes[i] for i in idx if mask[i])})
    return pd.DataFrame(rows)


def core_section(out: C.Output, rng: np.random.Generator, variant: str = "core") -> list[str]:
    """A fit whose prior holds no NDD set: which of its programmes carry the NDD lists?"""
    sub, prefix, title = VARIANTS[variant]
    fits = {(ds, chem): load(C.ns(ds, chem), out, sub) for ds, chem in C.STRATA}
    if all(v is None for v in fits.values()):
        return [f"**{title}** (no NDD set in the prior, to test the lists against programmes they did not shape): "
                f"not run yet -- slurm_spectra.sh with SPECTRA_PRIOR={variant}."]
    overview, support, enrich, pairs, ages, modnd = [], [], [], [], [], []
    mod_go = {}
    for ds in C.DATASETS:
        sets, steer = ndd_sets_for(ds, out), steering_sets(ds, out)
        if variant == "string":
            mods, go = string_module_sets(ds, sub, out)
            steer.update(mods)
            mod_go.update({(ds, m): t for m, t in go.items()})
            lc_v2, _ = C.cluster_expression(C.ns(ds, "v2"))
            md = module_ndd_content(mods, sets, lc_v2)
            if len(md):
                modnd.append(md.assign(dataset=ds))
        for chem in C.CHEMISTRIES:
            fit = fits.get((ds, chem))
            if fit is None:
                continue
            info, w, scores, _ = fit
            overview.append(info.reset_index().assign(dataset=ds, chemistry=chem, free=free_factors(info).to_numpy()))
            other = "v3" if chem == "v2" else "v2"
            lc_other, _ = C.cluster_expression(C.ns(ds, other))
            support.append(held_out(info, w, {}, lc_other, rng).assign(dataset=ds, chemistry=chem))
            lc_own, _ = C.cluster_expression(C.ns(ds, chem))
            e = ndd_enrichment(info, w, sets, lc_own, steer)
            if len(e):
                enrich.append(e.assign(dataset=ds, chemistry=chem))
            a = age_trends(scores, ds)
            if len(a):
                ages.append(a.assign(dataset=ds, chemistry=chem))
        if fits.get((ds, "v2")) is not None and fits.get((ds, "v3")) is not None:
            pairs.append(labelled_pairs(fits[(ds, "v2")], fits[(ds, "v3")], ds))
    overview = pd.concat(overview, ignore_index=True)
    support = pd.concat(support, ignore_index=True)
    support["q"] = C.bh(support.perm_p)
    pairs = pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame()
    enrich = pd.concat(enrich, ignore_index=True) if enrich else pd.DataFrame()
    ages = pd.concat(ages, ignore_index=True) if ages else pd.DataFrame()
    erep, acomb = pd.DataFrame(), pd.DataFrame()
    if len(enrich):
        enrich["q"] = np.nan
        for _, ix in enrich.groupby("dataset").groups.items():
            enrich.loc[ix, "q"] = C.bh(enrich.loc[ix, "p"])
    if len(enrich) and len(pairs):
        rows = []
        for p in pairs.itertuples():
            e2 = enrich[(enrich.dataset == p.dataset) & (enrich.chemistry == "v2") & (enrich.factor == p.factor_v2)]
            e3 = enrich[(enrich.dataset == p.dataset) & (enrich.chemistry == "v3") & (enrich.factor == p.factor_v3)]
            for r in e2.merge(e3, on="gene_set", suffixes=("_v2", "_v3")).itertuples():
                rows.append({"dataset": p.dataset, "pair": f"{p.factor_v2}~{p.factor_v3}", "scope": p.scope_v2,
                             "label_v2": p.label_v2, "label_v3": p.label_v3, "free": bool(r.free_v2 and r.free_v3),
                             "cosine": p.cosine, "gene_set": r.gene_set, "in_top_v2": r.in_top_v2,
                             "in_top_v3": r.in_top_v3, "fold_v2": r.fold_v2, "fold_v3": r.fold_v3,
                             "q_v2": r.q_v2, "q_v3": r.q_v3, "replicated": bool(r.q_v2 < 0.05 and r.q_v3 < 0.05),
                             "genes_in_both": "|".join(sorted(set(r.genes_v2.split("|")) & set(r.genes_v3.split("|"))
                                                              - {""}))})
        erep = pd.DataFrame(rows)
        acomb = combine_ages(pairs, ages) if len(ages) else pd.DataFrame()
    out.write(overview, f"{prefix}_factors", f"{title}: per stratum x factor (scope, label, free, top genes)")
    out.write(pairs, f"{prefix}_factor_pairs_v2_v3", f"{title}: mutual-best v2 / v3 factor pairs")
    out.write(support, f"{prefix}_held_out_support", f"{title}: top-gene coherence in the other chemistry's clusters")
    out.write(enrich, f"{prefix}_ndd_enrichment", f"{title}: per factor x NDD set, set genes among the top {N_TOP} "
              "vs the expression-matched expectation; exact p, BH")
    out.write(erep, f"{prefix}_ndd_enrichment_replicated", f"{title}: per v2 / v3 pair x NDD set")
    out.write(acomb, f"{prefix}_factor_age_trends_combined", f"{title}: per pair x class x region (or pooled): "
              "age trend, v2 x v3 combined")

    f = []
    ok = support[(support.q < 0.05) & (support.effect_vs_null_sd > 0)]
    what = {"core": "GO processes, cell-cycle panels and class markers only, more free factors",
            "string": "the core prior plus STRING network modules (experimental / database evidence)"}[variant]
    f.append(f"**{title}** ({what}; the NDD lists are tested against programmes they did not shape): " + "; ".join(
                 f"{ds} {ch}: {len(g)} factors ({int(g.free.sum())} free), "
                 f"{len(ok[(ok.dataset == ds) & (ok.chemistry == ch)])} supported in the other chemistry's donors"
                 for (ds, ch), g in overview.groupby(["dataset", "chemistry"])) + "; v2 / v3 pairs: "
             + ", ".join(f"{ds} {len(g)}" for ds, g in pairs.groupby("dataset")) + ".")
    if len(erep):
        rp = erep[erep.replicated].copy()
        lead = {(r.dataset, r.factor): ", ".join(str(r.top_genes).split("|")[:5])
                for r in overview[overview.chemistry == "v2"].itertuples()}
        sup = {(r.dataset, r.chemistry, r.factor): r.effect_vs_null_sd for r in support.itertuples()}
        if len(rp):
            rp["fmin"] = rp[["fold_v2", "fold_v3"]].min(axis=1)
            parts = []
            for (ds, pair), g in rp.sort_values("fmin", ascending=False).groupby(["dataset", "pair"], sort=False):
                r0 = g.iloc[0]
                f2_, f3_ = pair.split("~")
                kind = "free" if r0.free else r0.label_v2
                trend = ""
                if len(acomb):
                    t = acomb[(acomb.dataset == ds) & (acomb.pair == pair) & (acomb.tier.fillna("") != "")]
                    if len(t):
                        t = t.sort_values("combined_p").iloc[0]
                        trend = f"; activity {t.direction} in {t.region} {t.cell_class} ({t.tier})"
                parts.append(f"{ds} {pair} ({r0.scope}, {kind}; top {lead.get((ds, f2_), '')}; held-out "
                             f"{sup.get((ds, 'v2', f2_), np.nan):+.0f}/{sup.get((ds, 'v3', f3_), np.nan):+.0f} SD{trend}): "
                             + ", ".join(f"{x.gene_set} x{x.fold_v2:.1f}/x{x.fold_v3:.1f}" for x in g.itertuples())
                             + f" [{', '.join(str(g.iloc[0].genes_in_both).split('|')[:10])}]")
            f.append(f"**Programmes carrying NDD genes, in both chemistries** ({title.lower()}; set genes among the top "
                     f"{N_TOP} vs the expression-matched expectation -- for a steered programme, genes of its own prior "
                     "set left out -- q < 0.05 in v2 and v3; fold v2/v3; shared NDD genes): " + "; ".join(parts[:MAX_LISTED * 2]) + ".")
        else:
            f.append(f"**Programmes carrying NDD genes, in both chemistries** ({title.lower()}): none.")
    if modnd:
        modnd = pd.concat(modnd, ignore_index=True)
        modnd["q"] = np.nan
        for _, ix in modnd.groupby("dataset").groups.items():
            modnd.loc[ix, "q"] = C.bh(modnd.loc[ix, "p"])
        # a module's programme replicates when both fits hold a factor labelled with it, paired v2 / v3
        rep_mods = {}
        for p in pairs.itertuples():
            if str(p.label_v2).startswith("string:") and p.label_v2 == p.label_v3:
                s2 = support[(support.dataset == p.dataset) & (support.chemistry == "v2") & (support.factor == p.factor_v2)]
                s3 = support[(support.dataset == p.dataset) & (support.chemistry == "v3") & (support.factor == p.factor_v3)]
                ok2 = bool(len(s2) and s2.q.iloc[0] < 0.05 and s2.effect_vs_null_sd.iloc[0] > 0)
                ok3 = bool(len(s3) and s3.q.iloc[0] < 0.05 and s3.effect_vs_null_sd.iloc[0] > 0)
                rep_mods[(p.dataset, p.label_v2)] = (f"{p.factor_v2}~{p.factor_v3}", ok2 and ok3)
        modnd["programme_pair"] = [rep_mods.get((r.dataset, r.module), ("", False))[0] for r in modnd.itertuples()]
        modnd["programme_supported"] = [rep_mods.get((r.dataset, r.module), ("", False))[1] for r in modnd.itertuples()]
        modnd["best_go_term"] = [mod_go.get((r.dataset, r.module), "") for r in modnd.itertuples()]
        out.write(modnd, f"{prefix}_module_ndd_content", "Per STRING module x NDD set: set genes in the module vs the "
                  "expression-matched expectation; whether the module's programme replicates v2 / v3 and is supported "
                  "in held-out donors")
        n_rep = sum(1 for v in rep_mods.values() if v[1])
        hit = modnd[(modnd.q < 0.05) & modnd.programme_supported].sort_values("p")
        f.append(f"**STRING modules rich in NDD genes that form a replicated programme** ({len(rep_mods)} module "
                 f"programmes paired v2 / v3, {n_rep} also supported in held-out donors; NDD genes in the module vs the "
                 "expression-matched expectation, q < 0.05): " + ("; ".join(
                     f"{r.dataset} {r.module} ({r.best_go_term}; {r.programme_pair}) {r.gene_set} x{r.fold:.1f} "
                     f"[{', '.join(str(r.genes).split('|')[:8])}]" for r in hit.head(MAX_LISTED * 2).itertuples())
                     or "none") + ".")
    return f


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
    overview, pairs_all, support, ages, refine = [], [], [], [], []
    tops = {}
    for ds in C.DATASETS:
        f2, f3 = fits.get((ds, "v2")), fits.get((ds, "v3"))
        for chem, fit in (("v2", f2), ("v3", f3)):
            if fit is None:
                continue
            info, w, scores, pri = fit
            overview.append(info.reset_index().assign(dataset=ds, chemistry=chem, free=free_factors(info).to_numpy()))
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
        if f2 is not None and f3 is not None:
            p = labelled_pairs(f2, f3, ds)
            pairs_all.append(p)
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
    acomb = combine_ages(pairs, ages) if len(pairs) and len(ages) else pd.DataFrame()
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
    f.extend(core_section(out, rng, "core"))
    f.extend(core_section(out, rng, "string"))
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
            within = acomb[acomb.region != "all regions"]
            status = []
            for r in pooled.itertuples():
                w_ = within[(within.dataset == r.dataset) & (within.pair == r.pair) & (within.cell_class == r.cell_class)]
                sign = np.sign(r.stouffer_z)
                if not len(w_):
                    status.append("not testable within a region")
                elif ((w_.tier.fillna("") != "") & (np.sign(w_.stouffer_z) == sign)).any():
                    status.append("also within a region")
                elif ((np.sign(w_.effect_v2) == -sign) & (np.sign(w_.effect_v3) == -sign)).any():
                    status.append("reversed within " + ", ".join(
                        w_[(np.sign(w_.effect_v2) == -sign) & (np.sign(w_.effect_v3) == -sign)].region))
                elif ((np.sign(w_.effect_v2) == sign) & (np.sign(w_.effect_v3) == sign)).any():
                    status.append("same direction within " + ", ".join(
                        w_[(np.sign(w_.effect_v2) == sign) & (np.sign(w_.effect_v3) == sign)].region) + ", weaker")
                else:
                    status.append("chemistries disagree within regions")
            pooled = pooled.assign(within_region=status)
            acomb.loc[pooled.index, "within_region"] = status
            counts = pooled.within_region.str.extract(r"^(also|reversed|same direction|not testable|chemistries)")[0]
            f.append(f"**Pooled over regions**: {len(pooled)} replicated or supported trends -- "
                     + ", ".join(f"{n} {lab}" for lab, n in counts.value_counts().items())
                     + ". Reversed within a region (the pooled trend likely follows which regions were sampled at "
                     "each age): " + ("; ".join(
                         f"{r.dataset} {r.label} [{r.pair}] in {r.cell_class} {r.direction} ({r.within_region})"
                         for r in pooled[pooled.within_region.str.startswith("reversed")].head(MAX_LISTED).itertuples())
                         or "none") + ".")
    out.write(acomb, "factor_age_trends_combined", "Per factor pair x class x region (or pooled): v2 x v3 "
              "combined; tier; for pooled trends, how they behave within regions")
    out.summary(
        TITLE,
        "Which Spectra programmes are supported by donors the fit never saw, which replicate between the v2 and v3 "
        "fits, what did Spectra add to or drop from the NDD lists, and how does programme activity change with age?",
        ["Fits: scripts 24-25 (stratified ~25k-cell subsample per stratum; prior = gene lists, seed NDD and cell-cycle "
         "panels, core GO processes, class markers; extra free factors).",
         f"Replication: mutual-best cosine of gene weights >= {MATCH_MIN}. Held-out support: mean pairwise Spearman of a "
         f"factor's top {N_TOP} genes across the other chemistry's cluster pseudobulks vs {N_RANDOM} random sets matched "
         "on level x spread; added genes' correlation with kept prior genes vs matched random genes; BH.",
         f"NDD enrichment (core-prior fit, no NDD set in its prior): NDD set genes among each factor's top {N_TOP} vs "
         "the expectation from each gene's level x spread bin, exact Poisson-binomial p, a steered factor's own "
         "prior-set genes left out; BH per dataset; replicated when q < 0.05 for both factors of a v2 / v3 pair.",
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
        ["Compare with scHPF or cNMF (programmes fitted without any prior)."])


if __name__ == "__main__":
    main()
