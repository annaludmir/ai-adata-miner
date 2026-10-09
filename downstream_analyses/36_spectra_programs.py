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
     combined across chemistries (signed Stouffer, BH, tiered).
  E. NDD lists: for list-labelled factor pairs, genes added in both fits and
     list genes dropped in both (in the model but outside the top 200).

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
N_KEEP = 200
MATCH_MIN = 0.5
N_RANDOM = 500
MIN_AGES = 5
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
    s = scores[scores.grouping == "cell_class_x_age"].copy()
    parts = s.group.str.split(r" \| ", regex=True)
    s["cell_class"], s["age"] = parts.str[0], pd.to_numeric(parts.str[1], errors="coerce")
    s = s[(s.n_cells >= 20) & s.age.notna()]
    s = s[~s.age.map(lambda a: C.excluded(ds, age=a))]
    rows = []
    for cls, g in s.groupby("cell_class"):
        M = g.pivot_table(index="factor", columns="age", values="mean_score").dropna(axis=1)
        if M.shape[1] < MIN_AGES:
            continue
        x = M.columns.to_numpy(float)
        rho = C.spearman_rows(M.to_numpy(float), x)
        p, _ = C.spearman_perm_p(rho, x)
        rows.append(pd.DataFrame({"factor": M.index, "cell_class": cls, "n_ages": M.shape[1],
                                  "rho_vs_age": rho, "perm_p": p}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


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
                    k2 = set(f2[1].loc[r.factor_v2].sort_values(ascending=False).index[:N_KEEP])
                    k3 = set(f3[1].loc[r.factor_v3].sort_values(ascending=False).index[:N_KEEP])
                    members = set(C.mapped_lists(ds, collapse=False).get(r.label_v2.split(":", 1)[1], []))
                    in_model = members & set(f2[1].columns) & set(f3[1].columns)
                    refine.append({"dataset": ds, "gene_list": r.label_v2, "factor_v2": r.factor_v2,
                                   "factor_v3": r.factor_v3, "cosine": r.cosine,
                                   "list_genes_in_model": len(in_model),
                                   "added_in_both": "|".join(sorted((t2 & t3) - members)),
                                   "dropped_in_both": "|".join(sorted(in_model - k2 - k3)),
                                   "n_added_in_both": len((t2 & t3) - members),
                                   "n_dropped_in_both": len(in_model - k2 - k3)})
    overview = pd.concat(overview, ignore_index=True)
    pairs = pd.concat(pairs_all, ignore_index=True) if pairs_all else pd.DataFrame()
    support = pd.concat(support, ignore_index=True)
    ages = pd.concat(ages, ignore_index=True) if ages else pd.DataFrame()
    refine = pd.DataFrame(refine)
    acomb = pd.DataFrame()
    if len(pairs) and len(ages):
        rows = []
        for p in pairs.itertuples():
            for chem, f in (("v2", p.factor_v2), ("v3", p.factor_v3)):
                a = ages[(ages.dataset == p.dataset) & (ages.chemistry == chem) & (ages.factor == f)]
                for r in a.itertuples():
                    rows.append({"dataset": p.dataset, "pair": f"{p.factor_v2}~{p.factor_v3}",
                                 "label": p.label_v2 if p.label_v2 == p.label_v3 else f"{p.label_v2} / {p.label_v3}",
                                 "cell_class": r.cell_class, "chemistry": chem, "rho_vs_age": r.rho_vs_age,
                                 "perm_p": r.perm_p, "n_ages": r.n_ages})
        if rows:
            acomb = C.combine_chemistries(pd.DataFrame(rows), ["dataset", "pair", "label", "cell_class"],
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
    out.write(refine, "list_refinement", "Per list-labelled factor pair: genes added in both fits, list genes dropped in both")

    f = []
    trained = {}
    for ds, chem in C.STRATA:
        fs = C.EXPORTS / C.ns(ds, chem) / "25_spectra" / "fit_summary.csv"
        if fs.exists():
            out.used(f"{C.ns(ds, chem)}/25_spectra/fit_summary.csv")
            r = pd.read_csv(fs).iloc[0]
            if pd.notna(r.epochs_run):
                trained[(ds, chem)] = (f", {int(r.epochs_run)} of {int(r.epochs_max)} epochs"
                                       + (" (STOPPED EARLY)" if r.epochs_run < r.epochs_max else ""))
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
                     + ("; ".join(f"{r.dataset} {r.chemistry} {r.label}: {r.n_added} added, r {r.added_vs_kept_r:.2f} vs "
                                  f"{r.added_null_mean:.2f} ({', '.join(str(r.added_genes).split('|')[:8])})"
                                  for r in good.head(MAX_LISTED).itertuples()) if len(good) else "none") + ".")
    ref = refine[(refine.n_added_in_both > 0) | (refine.n_dropped_in_both > 0)] if len(refine) else refine
    if len(ref):
        f.append("**List refinement consistent across both fits**: " + "; ".join(
            f"{r.dataset} {r.gene_list}: +{r.n_added_in_both} added ({', '.join(r.added_in_both.split('|')[:8])}), "
            f"-{r.n_dropped_in_both} of {r.list_genes_in_model} dropped" for r in ref.itertuples()) + ".")
    if "best_08_module_overlap" in overview:
        nm = overview[(overview.label == "new") & (overview.best_08_module_overlap < 0.2)]
        f.append(f"**New programmes vs the 08 co-expression modules**: {int((overview.label == 'new').sum())} new "
                 f"factors, {len(nm)} overlapping no robust 08 module (top-gene overlap < 0.2)"
                 + (": " + "; ".join(f"{r.dataset} {r.chemistry} {r.factor} "
                                     f"({', '.join(str(r.top_genes).split('|')[:6])})"
                                     for r in nm.head(MAX_LISTED).itertuples()) if len(nm) else "") + ".")
    if len(acomb):
        rep = acomb[acomb.tier != ""]
        f.append("**Programme activity changing with age** (replicated pairs; rho v2/v3): " + ("; ".join(
            f"{r.dataset} {r.label} in {r.cell_class} {r.direction} ({r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f}; {r.tier})"
            for r in rep.head(MAX_LISTED * 2).itertuples()) if len(rep) else "none replicated") + ".")
    out.summary(
        TITLE,
        "Which Spectra programmes are supported by donors the fit never saw, which replicate between the v2 and v3 "
        "fits, what did Spectra add to or drop from the NDD lists, and how does programme activity change with age?",
        ["Fits: scripts 24-25 (stratified ~25k-cell subsample per stratum; prior = gene lists, seed NDD and cell-cycle "
         "panels, core GO processes, class markers; extra free factors).",
         f"Replication: mutual-best cosine of gene weights >= {MATCH_MIN}. Held-out support: mean pairwise Spearman of a "
         f"factor's top {N_TOP} genes across the other chemistry's cluster pseudobulks vs {N_RANDOM} random sets matched "
         "on level x spread; added genes' correlation with kept prior genes vs matched random genes; BH.",
         f"Age: Spearman of mean cell score per class x age (>= {MIN_AGES} ages); replicated pairs combined (signed "
         "Stouffer, BH, tiered). Lists: added = top-50 in both fits and not in the list; dropped = in the model but "
         f"outside the top {N_KEEP} in both."],
        f,
        ["A prior-steered factor can echo its prior whatever the data; only the held-out tests and the v2 / v3 "
         "replication speak to support.",
         "Fits use a subsample (~25k cells per stratum), so rare cell types contribute few cells.",
         "The gpu backend (Spectra_gpu, minibatched) is marked by its authors as in development; the cpu backend is "
         "the published model."],
        ["Compare with scHPF or cNMF (data-driven programmes) and with STRING-expanded priors."])


if __name__ == "__main__":
    main()
