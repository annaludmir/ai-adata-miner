#!/usr/bin/env python3
"""08 - Data-driven co-expression modules: which genes work as groups, unprompted?

Question: without starting from any list, which groups of genes rise and fall
together across the fine clusters of each file, robustly enough to be found
again in an independent donor set? Where does each module peak, what are its
hub genes, does it change with age inside cell classes -- and which gene lists
(user lists, seed NDD panels) concentrate in which module?

Method, per dataset and co-expression context (as in 07: across_clusters for
identity programmes, within_class for co-variation inside cell types):
  * Genes expressed in both chemistries; the 3,000 most variable in the
    context (mean of the v2 and v3 variance ranks).
  * Spearman correlation per chemistry; consensus = Fisher-z mean of the two.
  * Modules: average-linkage clustering of 1 - r, cut at the 99th percentile
    of the consensus correlations (kept within 0.3-0.5: within-class
    residuals correlate far less than raw profiles); any module above 250
    genes is re-clustered on its own at a stricter cut (+0.1, up to 0.8) -- a
    single cut leaves two giant neuron-vs-progenitor halves. Modules under 20
    genes are left unassigned.
  * Robustness: the same procedure runs on v2 alone and v3 alone. A consensus
    module is robust when its best match among the v2-only modules holds up
    in v3 and its best match among the v3-only modules holds up in v2 (held-out
    coherence vs matched random gene sets, BH), each match with Jaccard >= 0.3.
    This is the standard reference/test module-preservation design, with
    disjoint donors as reference and test.
  * Eigengene = mean z-scored member profile per cluster; hub genes by
    correlation with it (kME, minimum over chemistries).
  * Age dynamics: module score per (cell class, age point) from TMM log CPM,
    Spearman with age (exact permutation), v2/v3 signed Stouffer, tiered.
  * List enrichment: overlap of each user list / seed NDD panel with each
    module vs random sets matched on expression level x spread within the
    3,000-gene universe.

Inputs:
  csv_exports/<ds>__<chem>/09_pseudobulk/<finest clustering>__{pseudobulk_counts,group_summary}.csv
  csv_exports/<ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,detection_fraction,group_summary}.csv
  csv_exports/<ds>__<chem>/04_clusters/cluster_profile_<clustering>.csv
  csv_exports/<ds>__<chem>/11_panels/panel_coverage.csv
  gene lists: config.GENE_LISTS_DIR (optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

import _common as C

SLUG = "08_coexpression_modules"
N_VARIABLE = 3000
R_CUT_QUANTILE = 0.99   # cut = this quantile of the consensus gene-gene correlations ...
R_CUT_RANGE = (0.3, 0.5)  # ... kept within this range
R_STEP = 0.1
R_MAX = 0.8
MAX_MODULE = 250
MIN_MODULE = 20
MIN_JACCARD = 0.3
N_RANDOM = 1000
N_HUBS = 10
CONTEXTS = ("across_clusters", "within_class")
ANNOTATION_GROUPS = ("cell_cycle", "marker", "patterning")
PANEL_LABELS = {"cell_cycle:g2m_phase": "G2/M phase", "cell_cycle:s_phase": "S phase"}


def pretty_panel(name: str) -> str:
    """'cell_cycle:g2m_phase' -> 'G2/M phase'; 'marker:radial_glia' -> 'radial glia markers'."""
    if name in PANEL_LABELS:
        return PANEL_LABELS[name]
    grp, panel = name.split(":", 1)
    suffix = {"marker": "markers", "patterning": "patterning"}.get(grp, "")
    return f"{panel.replace('_', ' ')} {suffix}".strip()
MAX_LISTED = 8


# ---------------------------------------------------------------------------
# module detection
# ---------------------------------------------------------------------------
def detect_modules(Cm: np.ndarray, genes: np.ndarray, r_cut: float,
                   prefix: str = "") -> dict[str, list[str]]:
    """{path label: genes} -- average linkage at r_cut, oversized modules re-split."""
    if len(genes) < MIN_MODULE:
        return {}
    D = np.clip(1 - Cm, 0, 2)
    np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), t=1 - r_cut,
                   criterion="distance")
    sizes = pd.Series(lab).value_counts()
    out = {}
    for i, (c, n) in enumerate(sizes[sizes >= MIN_MODULE].items(), start=1):
        sel = lab == c
        label = f"{prefix}{i}"
        if n > MAX_MODULE and r_cut + R_STEP <= R_MAX + 1e-9:
            sub = detect_modules(Cm[np.ix_(sel, sel)], genes[sel], r_cut + R_STEP, label + ".")
            if len(sub) >= 2:
                out.update(sub)
                continue
        out[label] = list(genes[sel])
    return out


def consensus(R2: np.ndarray, R3: np.ndarray) -> np.ndarray:
    Cm = np.tanh((np.arctanh(np.clip(R2, -0.999, 0.999)) + np.arctanh(np.clip(R3, -0.999, 0.999))) / 2)
    np.fill_diagonal(Cm, 1.0)
    return Cm


def held_out_p(Zn: np.ndarray, bins: np.ndarray, idx: np.ndarray, rng) -> tuple[float, float, float]:
    """(coherence, null mean, one-sided p) of a gene set in a matrix it was not found in."""
    k = idx.size
    obs = float(C.coherence_from_sums(Zn[idx].sum(axis=0)[None, :], k)[0])
    null = C.coherence_from_sums(C.set_sums(Zn, C.matched_sets(bins, idx, N_RANDOM, rng)), k)
    return obs, float(null.mean()), float(((null >= obs).sum() + 1) / (N_RANDOM + 1))


def matched_p(null: np.ndarray, obs: float) -> float:
    """One-sided p of an overlap against its matched-null draws.

    Empirical while any draw reaches the observed value; beyond the null's
    resolution (no draw reaches it), a normal approximation from the null's
    mean and SD. Otherwise every strong overlap ties at 1/(N+1), and with a few
    hundred tests BH cannot pass even a complete overlap.
    """
    k = int((null >= obs).sum())
    if k > 0:
        return float((k + 1) / (null.size + 1))
    sd = float(null.std())
    if sd <= 0:
        return float(1 / (null.size + 1))
    return float(max(stats.norm.sf((obs - null.mean()) / sd), 1e-300))


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


# ---------------------------------------------------------------------------
def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)

    data = {}
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        grouping = C.COEXPR_GROUPING[ds]
        out.used(f"{n}/09_pseudobulk/{grouping}__pseudobulk_counts.csv",
                 f"{n}/04_clusters/cluster_profile_{grouping.replace('cluster_', '', 1)}.csv",
                 f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv",
                 f"{n}/11_panels/panel_coverage.csv")
        lc, annot = C.cluster_expression(n)
        data[(ds, chem)] = (lc, annot, C.coexpr_contexts(lc, annot))

    # gene sets to test against modules: user lists + seed NDD panels
    gene_sets = {}
    for ds in C.DATASETS:
        sets = {f"list:{k}": v for k, v in C.mapped_lists(ds).items()} if C.gene_lists() else {}
        pan = C.panels(C.ns(ds, "v2"))
        for pname, g in pan[pan.panel_group == "ndd"].groupby("panel"):
            sets[f"seed:{pname}"] = list(g.gene)
        gene_sets[ds] = sets
    # seed reference panels used only to label modules (cell-cycle phase,
    # cell-class markers, patterning) -- not reported as findings themselves
    annot_sets = {}
    for ds in C.DATASETS:
        pan = C.panels(C.ns(ds, "v2"))
        annot_sets[ds] = {f"{grp}:{pname}": list(g.gene)
                          for (grp, pname), g in pan[pan.panel_group.isin(ANNOTATION_GROUPS)]
                          .groupby(["panel_group", "panel"])}

    mod_rows, mem_rows, age_rows, enr_rows, match_rows, ann_rows = [], [], [], [], [], []
    eig_tables = {}
    for ds in C.DATASETS:
        lcs = {c: data[(ds, c)][0] for c in C.CHEMISTRIES}
        shared = lcs["v2"].index.intersection(lcs["v3"].index)
        pos = {c: {g: i for i, g in enumerate(lcs[c].index)} for c in C.CHEMISTRIES}
        ages_tab = {c: C.class_age_logcpm(C.ns(ds, c)) for c in C.CHEMISTRIES}
        for cname in CONTEXTS:
            if not all(cname in data[(ds, c)][2] for c in C.CHEMISTRIES):
                continue
            ctx = {c: data[(ds, c)][2][cname] for c in C.CHEMISTRIES}
            # variable genes in this context: mean variance rank over chemistries
            rk = []
            for c in C.CHEMISTRIES:
                V = C.context_values(lcs[c], data[(ds, c)][1], cname)
                rk.append(V.var(axis=1).loc[shared].rank(ascending=False))
            universe = pd.concat(rk, axis=1).mean(axis=1).sort_values().index[:N_VARIABLE].to_numpy()
            rows = {c: np.array([pos[c][g] for g in universe]) for c in C.CHEMISTRIES}
            R = {c: ctx[c][0][rows[c]] @ ctx[c][0][rows[c]].T for c in C.CHEMISTRIES}
            Cm = consensus(R["v2"], R["v3"])
            # Correlation strength differs by dataset and context (within-class
            # residuals correlate far less), so the cut follows the data: the
            # strongest ~1% of gene pairs, bounded to a sensible range.
            r_cut = float(np.clip(np.quantile(Cm[np.triu_indices_from(Cm, 1)], R_CUT_QUANTILE),
                                  *R_CUT_RANGE))
            mods = detect_modules(Cm, universe, r_cut)
            single = {c: detect_modules(R[c], universe, r_cut) for c in C.CHEMISTRIES}
            C.log(f"  {ds} / {cname}: cut r >= {r_cut:.2f}; {len(mods)} consensus modules covering "
                  f"{sum(len(v) for v in mods.values())} of {len(universe)} genes; "
                  + ", ".join(f"{len(single[c])} in {c} alone" for c in C.CHEMISTRIES))

            # held-out preservation of single-chemistry modules
            pres = {}
            for disc, held in (("v2", "v3"), ("v3", "v2")):
                for lab, genes in single[disc].items():
                    idx = np.array([pos[held][g] for g in genes])
                    obs, nm, p = held_out_p(ctx[held][0], ctx[held][1], idx, rng)
                    pres[(disc, lab)] = (obs, nm, p)
            keys = list(pres)
            qv = dict(zip(keys, C.bh([pres[k][2] for k in keys]))) if keys else {}

            # name consensus modules by size, annotate
            order = sorted(mods, key=lambda k: -len(mods[k]))
            prefix = {"cortex": "CX", "human_dev": "HD"}[ds] + ("" if cname == "across_clusters" else "w")
            names = {lab: f"{prefix}{i:02d}" for i, lab in enumerate(order, start=1)}
            scores = {}
            for lab in order:
                genes = mods[lab]
                name = names[lab]
                row = {"dataset": ds, "context": cname, "module": name, "tree_path": lab,
                       "r_cut": r_cut, "n_genes": len(genes)}
                robust = True
                for disc, held in (("v2", "v3"), ("v3", "v2")):
                    best = max(single[disc], key=lambda k: jaccard(genes, single[disc][k]), default=None)
                    j = jaccard(genes, single[disc][best]) if best else 0.0
                    obs, nm, p = pres.get((disc, best), (np.nan, np.nan, np.nan))
                    q = qv.get((disc, best), np.nan)
                    row.update({f"match_{disc}_only_jaccard": j,
                                f"match_{disc}_only_heldout_coherence": obs,
                                f"match_{disc}_only_heldout_null": nm,
                                f"match_{disc}_only_heldout_q": q})
                    robust &= bool(j >= MIN_JACCARD and np.isfinite(q) and q < 0.05)
                row["robust"] = robust
                # eigengene, hubs, annotation per chemistry
                kme = {}
                for c in C.CHEMISTRIES:
                    Zn, _, cols, Zz = ctx[c]
                    idx = np.array([pos[c][g] for g in genes])
                    e = Zz[idx].mean(axis=0)
                    scores[(name, c)] = pd.Series(e, index=cols)
                    ez = C.unit_rank_rows(e[None, :])[0]
                    kme[c] = pd.Series(Zn[idx] @ ez, index=genes)
                    a = data[(ds, c)][1].loc[cols]
                    # Where the genes are expressed comes from their raw profiles: in the
                    # within-class context every class mean is zero by construction.
                    _, _, cols_raw, Zz_raw = data[(ds, c)][2]["across_clusters"]
                    a_raw = data[(ds, c)][1].loc[cols_raw]
                    e_raw = pd.Series(Zz_raw[idx].mean(axis=0), index=cols_raw)
                    by_cls = e_raw.groupby(a_raw["dominant_cell_class"].to_numpy()).mean()
                    by_cls = by_cls[a_raw["dominant_cell_class"].value_counts().reindex(by_cls.index) >= 3]
                    row[f"peak_class_{c}"] = by_cls.idxmax() if len(by_cls) else ""
                    if ds == "human_dev":
                        by_reg = e_raw.groupby(a_raw["dominant_region"].to_numpy()).mean()
                        by_reg = by_reg[a_raw["dominant_region"].value_counts().reindex(by_reg.index) >= 3]
                        row[f"peak_region_{c}"] = by_reg.idxmax() if len(by_reg) else ""
                    top = pd.Series(e, index=cols).sort_values(ascending=False).head(3)
                    row[f"top_clusters_{c}"] = "; ".join(
                        f"{x} ({a.loc[x, 'dominant_cell_class']}, {a.loc[x, 'dominant_region']}, "
                        f"{a.loc[x, 'age_pcw_median']:g} pcw)" for x in top.index)
                    row[f"rho_with_cluster_age_{c}"] = float(
                        stats.spearmanr(e, a["age_pcw_median"].to_numpy(float)).statistic)
                kmin = pd.concat(kme, axis=1).min(axis=1).sort_values(ascending=False)
                row["hub_genes"] = "|".join(kmin.index[:N_HUBS])
                row["genes"] = "|".join(kmin.index)
                mod_rows.append(row)
                for g in genes:
                    mem_rows.append({"dataset": ds, "context": cname, "module": name, "gene": g,
                                     "kme_v2": kme["v2"][g], "kme_v3": kme["v3"][g],
                                     "kme_min": kmin[g], "hub": g in set(kmin.index[:N_HUBS])})

                # age dynamics within classes
                for c in C.CHEMISTRIES:
                    for cls, (lcc, _, ages, _) in ages_tab[c].items():
                        if len(ages) < 5:
                            continue
                        present = [g for g in genes if g in lcc.index]
                        if len(present) < MIN_MODULE // 2:
                            continue
                        sc = C.zscore_rows(lcc.loc[present].to_numpy(float)).mean(axis=0)
                        rho = C.spearman_rows(sc[None, :], ages)
                        p, exact = C.spearman_perm_p(rho, ages, X=sc[None, :])
                        age_rows.append({"dataset": ds, "context": cname, "module": name,
                                         "chemistry": c, "cell_class": cls, "n_ages": len(ages),
                                         "rho": float(rho[0]), "perm_p": float(p[0])})

            # list / panel enrichment against this context's modules
            uni = pd.Index(universe)
            lc2u = lcs["v2"].loc[universe]
            bins = C.expression_bins(lc2u)
            member = np.zeros((len(order), len(universe)), dtype=float)
            for i, lab in enumerate(order):
                member[i, uni.get_indexer(mods[lab])] = 1.0
            for sets, rows_out in ((gene_sets[ds], enr_rows), (annot_sets[ds], ann_rows)):
              for sname, sgenes in sets.items():
                idx = uni.get_indexer([g for g in dict.fromkeys(sgenes) if g in uni])
                idx = idx[idx >= 0]
                if idx.size < 5:
                    continue
                obs = member[:, idx].sum(axis=1)
                rand = C.matched_sets(bins, idx, N_RANDOM, rng)
                null = np.stack([member[:, r].sum(axis=1) for r in rand], axis=1)
                for i, lab in enumerate(order):
                    exp = float(null[i].mean())
                    rows_out.append({"dataset": ds, "context": cname, "gene_set": sname,
                                     "module": names[lab], "set_genes_in_universe": int(idx.size),
                                     "module_size": len(mods[lab]), "overlap": int(obs[i]),
                                     "expected_matched": exp,
                                     "fold": obs[i] / exp if exp > 0 else np.nan,
                                     "p_matched": matched_p(null[i], obs[i]),
                                     "genes": "|".join(uni[idx][member[i, idx] > 0])})
            eig_tables[(ds, cname)] = scores

    mods_df = pd.DataFrame(mod_rows)
    # label modules from the reference panels they concentrate
    ann = pd.DataFrame(ann_rows)
    if not ann.empty:
        ann["q"] = np.nan
        for _, ix in ann.groupby(["dataset", "context"]).groups.items():
            ann.loc[ix, "q"] = C.bh(ann.loc[ix, "p_matched"])
        ann = ann.sort_values(["dataset", "context", "p_matched"])
        hits = ann[(ann.q < 0.05) & (ann.fold >= 2)]
        labels = (hits.groupby(["dataset", "context", "module"])["gene_set"]
                  .apply(lambda g: "; ".join(pretty_panel(x) for x in list(g)[:2])))
        mods_df["annotation"] = [labels.get((r.dataset, r.context, r.module), "")
                                 for r in mods_df.itertuples()]
        out.write(ann, "module_annotation",
                  "Overlap of each module with seed reference panels (cell-cycle phase, cell-class "
                  "markers, patterning) vs matched random sets; used to label modules")
    else:
        mods_df["annotation"] = ""
    mem_df = pd.DataFrame(mem_rows)
    out.write(mods_df, "modules",
              "Per module: size, robustness (independent discovery + held-out preservation), "
              "peak class/region, top clusters, hub genes, members")
    out.write(mem_df, "module_membership", "Gene -> module, with kME per chemistry and hub flag")

    age = pd.DataFrame(age_rows)
    age_comb = []
    for (ds, cname, mod, cls), g in age.groupby(["dataset", "context", "module", "cell_class"], sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        Z, pc = C.signed_stouffer([np.array([a.rho]), np.array([b.rho])],
                                  [np.array([a.perm_p]), np.array([b.perm_p])],
                                  [np.sqrt(a.n_ages), np.sqrt(b.n_ages)])
        age_comb.append({"dataset": ds, "context": cname, "module": mod, "cell_class": cls,
                         "rho_v2": a.rho, "p_v2": a.perm_p, "rho_v3": b.rho, "p_v3": b.perm_p,
                         "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
    age_comb = pd.DataFrame(age_comb)
    if not age_comb.empty:
        age_comb["combined_q"] = np.nan
        for _, ix in age_comb.groupby(["dataset", "context"]).groups.items():
            age_comb.loc[ix, "combined_q"] = C.bh(age_comb.loc[ix, "combined_p"])
        age_comb["tier"] = C.replication_tier(age_comb.rho_v2, age_comb.p_v2, age_comb.rho_v3,
                                              age_comb.p_v3, age_comb.combined_q)
        age_comb["direction"] = np.where(age_comb.stouffer_z > 0, "rises", "falls")
        age_comb = age_comb.sort_values(["dataset", "context", "combined_p"])
    out.write(age, "module_age_trends_per_stratum",
              "Module score vs age within each cell class, per chemistry (exact permutation p)")
    out.write(age_comb, "module_age_trends_combined",
              "v2 x v3 combined module age trends per class; tier replicated / supported")

    enr = pd.DataFrame(enr_rows)
    if not enr.empty:
        enr["q"] = np.nan
        for _, ix in enr.groupby(["dataset", "context"]).groups.items():
            enr.loc[ix, "q"] = C.bh(enr.loc[ix, "p_matched"])
        enr = enr.sort_values(["dataset", "context", "p_matched"])
    out.write(enr, "module_gene_set_enrichment",
              "Overlap of each user list / seed NDD panel with each module vs matched random sets")

    # cross-dataset module matching (same donors: reproducibility of processing)
    for cname in CONTEXTS:
        a = mods_df[(mods_df.dataset == "cortex") & (mods_df.context == cname)]
        b = mods_df[(mods_df.dataset == "human_dev") & (mods_df.context == cname)]
        for r in a.itertuples():
            ga = set(r.genes.split("|"))
            if b.empty:
                continue
            js = b.genes.map(lambda x: jaccard(ga, x.split("|")))
            j = js.idxmax()
            match_rows.append({"context": cname, "cortex_module": r.module,
                               "human_dev_module": b.loc[j, "module"], "jaccard": js.max(),
                               "shared_genes": len(ga & set(b.loc[j, "genes"].split("|")))})
    matches = pd.DataFrame(match_rows)
    out.write(matches, "cortex_vs_human_dev_module_match",
              "Best-matching human_dev module for each cortex module (Jaccard; same donors)")

    figures(out, mods_df, eig_tables, data, enr)

    # ---- findings -------------------------------------------------------------
    f = []
    for (ds, cname), g in mods_df.groupby(["dataset", "context"], sort=False):
        rob = g[g.robust]
        items = []
        for r in rob.head(MAX_LISTED).itertuples():
            peak = r.peak_class_v2 if r.peak_class_v2 == r.peak_class_v3 else f"{r.peak_class_v2}/{r.peak_class_v3}"
            reg = ""
            if ds == "human_dev":
                pr2, pr3 = getattr(r, "peak_region_v2", ""), getattr(r, "peak_region_v3", "")
                reg = f", {pr2}" if pr2 == pr3 else f", {pr2}/{pr3}"
            trends = age_comb[(age_comb.dataset == ds) & (age_comb.context == cname)
                              & (age_comb.module == r.module) & (age_comb.tier == "replicated")]
            ttxt = ("; " + ", ".join(f"{t.direction} with age in {t.cell_class}"
                                     for t in trends.itertuples())) if len(trends) else ""
            sets = enr[(enr.dataset == ds) & (enr.context == cname) & (enr.module == r.module)
                       & (enr.q < 0.05) & (enr.fold >= 2)] if not enr.empty else enr
            stxt = ("; holds " + ", ".join(f"{s.gene_set} ({s.overlap}, {s.fold:.1f}x)"
                                           for s in sets.itertuples())) if len(sets) else ""
            tag = f" [{r.annotation}]" if getattr(r, "annotation", "") else ""
            items.append(f"{r.module}{tag} ({r.n_genes} genes; peak {peak}{reg}; hubs "
                         f"{', '.join(r.hub_genes.split('|')[:5])}{ttxt}{stxt})")
        if len(rob) > MAX_LISTED:
            items.append(f"(+{len(rob) - MAX_LISTED} more robust modules in modules.csv)")
        f.append(f"**{ds} / {cname.replace('_', ' ')}: {len(g)} modules, {len(rob)} robust** "
                 "(re-discovered in each chemistry alone and preserved in the other; cut r >= "
                 f"{g.r_cut.iloc[0]:.2f})" + "".join(f"\n  - {it}" for it in items))
    if not enr.empty:
        sig = enr[(enr.q < 0.05) & (enr.fold >= 2)]
        robust_names = set(zip(mods_df.loc[mods_df.robust, "dataset"], mods_df.loc[mods_df.robust, "context"],
                               mods_df.loc[mods_df.robust, "module"]))
        sig = sig[[k in robust_names for k in zip(sig.dataset, sig.context, sig.module)]]
        if len(sig):
            by_set = sig.groupby("gene_set").apply(
                lambda d: "; ".join(f"{r.dataset} {r.module} ({r.overlap} genes, {r.fold:.1f}x)"
                                    for r in d.sort_values("p_matched").head(4).itertuples()),
                include_groups=False)
            f.append("**Gene sets concentrated in robust modules** (q < 0.05, >= 2x matched "
                     "expectation): " + " | ".join(f"{k}: {v}" for k, v in by_set.items()) + ".")
        else:
            f.append("**No user list or seed panel concentrates in a robust module** beyond "
                     "matched expectation (q < 0.05, >= 2x).")
    if not matches.empty:
        good = matches[matches.jaccard >= MIN_JACCARD]
        f.append(f"**cortex and human_dev share {len(good)} of {len(matches)} cortex modules** "
                 f"(Jaccard >= {MIN_JACCARD} with a human_dev module; same donors, so this is "
                 "reproducibility of processing and regional pooling, not replication).")

    out.summary(
        "Co-expression modules across fine clusters",
        "Without starting from a list: which groups of genes co-vary across cell clusters "
        "robustly enough to be found again in independent donors, where do they peak, how do "
        "they change with age, and which gene lists concentrate in them?",
        [f"Clusters as in 07 ({C.COEXPR_GROUPING['cortex']}, {C.COEXPR_GROUPING['human_dev']}; "
         f">= {C.COEXPR_MIN_CELLS} cells; donor-dominated dropped; log2 TMM-CPM; Spearman).",
         f"Universe: the {N_VARIABLE:,} most variable genes expressed in both chemistries, per "
         "context (across clusters; within cell classes after removing class means).",
         "Modules: average linkage on 1 - r of the v2/v3 consensus (Fisher z), cut at the "
         f"{R_CUT_QUANTILE:.0%} quantile of the consensus correlations, bounded to "
         f"{R_CUT_RANGE[0]}-{R_CUT_RANGE[1]} (within-class residuals correlate far less than "
         f"profiles do); modules > {MAX_MODULE} genes re-split at +{R_STEP} (to {R_MAX}); < "
         f"{MIN_MODULE} genes unassigned. Named CX/HD + number (w = within-class context).",
         "Robust = the best-matching module found in v2 alone is preserved in v3, and vice "
         f"versa (held-out coherence vs {N_RANDOM:,} matched random sets, BH q < 0.05; "
         f"Jaccard >= {MIN_JACCARD}).",
         "Age trends: module score (mean z of members over a class's age points, TMM log CPM) "
         "vs age, exact permutation p, v2/v3 signed Stouffer, tiered as elsewhere.",
         "Enrichment: overlap vs random sets matched on level x spread within the universe "
         "(empirical p; beyond the null's resolution, a normal approximation from its mean and "
         "SD); BH per dataset x context.",
         "Labels in [brackets]: the seed reference panels (cell-cycle phase, cell-class markers, "
         "patterning) a module concentrates (q < 0.05, >= 2x matched expectation; top two)."],
        f,
        ["Modules depend on the cut and the universe: they are a summary of correlation "
         "structure, not discrete biological units. The tree path in modules.csv shows which "
         "modules came from re-splitting one larger module.",
         "Cluster pseudobulks are not independent samples; preservation is judged against "
         "matched random gene sets in held-out donors, not by absolute correlation.",
         "Gene-set nulls ignore correlation among set members; replication and fold >= 2 are "
         "the guards.",
         "human_dev modules pool brain regions; cortex and human_dev share donors."],
        ["Annotate robust modules with GO / TF-target enrichment (needs annotation files).",
         "Compare module eigengenes between regions once a region x cluster export exists."])


def figures(out: C.Output, mods_df: pd.DataFrame, eig_tables: dict, data: dict,
            enr: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or mods_df.empty:
        return
    for (ds, cname), scores in eig_tables.items():
        g = mods_df[(mods_df.dataset == ds) & (mods_df.context == cname)]
        if g.empty:
            continue
        mats = []
        for c in C.CHEMISTRIES:
            a = data[(ds, c)][1]
            lc = data[(ds, c)][0]
            _, _, cols_raw, Zz_raw = data[(ds, c)][2]["across_clusters"]
            pos = {gname: i for i, gname in enumerate(lc.index)}
            cls_raw = a.loc[cols_raw, "dominant_cell_class"].to_numpy()
            m = {}
            for mod, genes in zip(g.module, g.genes):
                idx = [pos[x] for x in genes.split("|") if x in pos]
                m[mod] = pd.Series(Zz_raw[idx].mean(axis=0), index=cols_raw).groupby(cls_raw).mean()
            mats.append(pd.DataFrame(m).T)
        classes = sorted(set(mats[0].columns) & set(mats[1].columns))
        M = (mats[0][classes] + mats[1][classes]) / 2
        fig, ax = plt.subplots(figsize=(0.6 * len(classes) + 3, 0.32 * len(M) + 1.8))
        im = ax.imshow(M.to_numpy(float), cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
        ax.set_xticks(range(len(classes)), classes, rotation=60, ha="right", fontsize=7)
        labels = [f"{m}{' *' if r else ''} ({n})" for m, r, n in zip(g.module, g.robust, g.n_genes)]
        ax.set_yticks(range(len(M)), labels, fontsize=6)
        ax.set_title(f"{ds} / {cname.replace('_', ' ')}: where module genes are expressed "
                     "(mean of v2, v3; * robust)", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.7, label="mean member z-score by class")
        fig.tight_layout()
        out.figure(fig, f"modules_{ds}_{cname}", f"Module eigengene by cell class, {ds} / {cname}")
        plt.close(fig)
    if enr is None or enr.empty:
        return
    for (ds, cname), g in enr.groupby(["dataset", "context"]):
        h = g.pivot(index="gene_set", columns="module", values="fold")
        sig = g.pivot(index="gene_set", columns="module", values="q") < 0.05
        h = np.log2(h.astype(float).clip(lower=0.125))
        fig, ax = plt.subplots(figsize=(0.32 * h.shape[1] + 3, 0.35 * h.shape[0] + 1.8))
        im = ax.imshow(h.to_numpy(float), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
        for (i, j) in zip(*np.nonzero(sig.to_numpy())):
            ax.text(j, i, "*", ha="center", va="center", fontsize=7)
        ax.set_xticks(range(h.shape[1]), h.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(h.shape[0]), h.index, fontsize=6)
        ax.set_title(f"{ds} / {cname.replace('_', ' ')}: gene sets in modules (log2 fold vs "
                     "matched; * q < 0.05)", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.7)
        fig.tight_layout()
        out.figure(fig, f"module_enrichment_{ds}_{cname}", f"Gene-set enrichment in modules, {ds} / {cname}")
        plt.close(fig)


if __name__ == "__main__":
    main()
