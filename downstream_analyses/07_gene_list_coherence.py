#!/usr/bin/env python3
"""07 - Do the genes in each list work as a group? Co-expression coherence.

Question: across the fine clusters of each file, do a list's genes rise and
fall together more than random genes of the same expression level and
spread? Which members carry that coherence, does the list split into tighter
sub-modules, where do those sub-modules peak -- and which genes outside the
list track them closely enough to be candidate members?

Two co-expression contexts per stratum:
  across_clusters -- correlation over all clusters: genes that switch on and
                     off together between cell types and states.
  within_class    -- each gene's mean per dominant cell class removed first:
                     co-variation *inside* cell types (sub-type, region, age),
                     beyond the shared "these are all neuronal genes" signal.

Method (per stratum and context):
  * Cluster pseudobulks with >= 100 cells, not donor-dominated, -> log2
    TMM-CPM; genes >= 10 CPM in >= 3 clusters. Spearman correlation.
  * Coherence = mean pairwise correlation among the list's genes. Null:
    1,000 random sets drawing, member by member, from the same bin of mean
    level x spread (5 x 5 quantile bins), without repeats. One-sided p.
  * Per-gene connectivity = mean correlation with the other members, against
    that gene's correlation with the random sets.
  * v2 and v3 (disjoint donors) combined by Stouffer; replicated = both
    nominal and combined q < 0.05.
  * Sub-modules, cross-validated: each chemistry splits the list on its own
    correlations (average linkage, cut at the best silhouette; no split if no
    cut reaches 0.1), and every sub-module is re-tested in the OTHER chemistry,
    whose donors played no part in finding it. Held-out coherence alone is not
    enough: in development a random list split into groups that also held up
    in the other donor set, because the transcriptome itself has reproducible
    structure. So each sub-module is compared with the sub-modules the same
    procedure carves out of 50 random lists matched to the real one (held-out
    coherence, comparable size). Sub-modules found in both directions with
    overlapping members (Jaccard >= 0.3) are reported with their shared core.
  * Candidates: non-members whose correlation with a validated sub-module's
    eigengene is >= 0.6 in both chemistries -- hypotheses, not members.

Inputs:
  gene lists: config.GENE_LISTS_DIR (AIM_GENE_LISTS overrides)
  csv_exports/<ds>__<chem>/09_pseudobulk/<finest clustering>__{pseudobulk_counts,group_summary}.csv
  csv_exports/<ds>__<chem>/04_clusters/cluster_profile_<clustering>.csv
  csv_exports/<ds>__<chem>/11_panels/panel_coverage.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

import _common as C

SLUG = "07_gene_list_coherence"
TITLE = "Gene-list co-expression: do list genes work as a group?"
N_RANDOM = 1000
MIN_GENES = 5
MIN_SUBMODULE = 5
MAX_SUBMODULES = 8
MIN_SILHOUETTE = 0.1
N_RANDOM_LISTS = 50        # random matched lists split the same way (selection-aware null)
MIN_NULL_SUBMODULES = 20   # comparable random sub-modules needed for a p-value
CAND_R = 0.6
N_CAND = 30
CONTEXTS = ("across_clusters", "within_class")
MAX_LISTED = 6
HEATMAP_MAX_GENES = 300


def coherence_test(Zn: np.ndarray, bins: np.ndarray, idx: np.ndarray,
                   rng: np.random.Generator):
    """(coherence, null array, per-gene connectivity, per-gene p)."""
    k = idx.size
    S = Zn[idx].sum(axis=0)
    coh = float(C.coherence_from_sums(S[None, :], k)[0])
    rand = C.matched_sets(bins, idx, N_RANDOM, rng)
    sums = C.set_sums(Zn, rand)
    null = C.coherence_from_sums(sums, k)
    conn = (C.dot(Zn[idx], S) - 1.0) / (k - 1)
    conn_null = C.dot(Zn[idx], sums.T) / k                   # gene vs random set
    p_gene = ((conn_null >= conn[:, None]).sum(axis=1) + 1) / (N_RANDOM + 1)
    return coh, null, conn, p_gene


def one_sided_combine(p2: float, p3: float) -> float:
    z = (stats.norm.isf(np.clip(p2, 1e-300, 1)) + stats.norm.isf(np.clip(p3, 1e-300, 1))) / np.sqrt(2)
    return float(stats.norm.sf(z))


def tier(p2, p3, q) -> str:
    if q < 0.05 and p2 < 0.05 and p3 < 0.05:
        return "replicated"
    return "supported" if q < 0.05 else ""


def silhouette_precomputed(D: np.ndarray, labels: np.ndarray) -> float:
    """Mean silhouette from a distance matrix (vectorised; singletons score 0)."""
    labels = np.asarray(labels)
    uniq, inv = np.unique(labels, return_inverse=True)
    if len(uniq) < 2:
        return np.nan
    onehot = np.zeros((len(labels), len(uniq)))
    onehot[np.arange(len(labels)), inv] = 1.0
    counts = onehot.sum(axis=0)
    sums = C.dot(D, onehot)
    own_n = counts[inv]
    a = np.divide(sums[np.arange(len(labels)), inv], own_n - 1,
                  out=np.zeros(len(labels)), where=own_n > 1)
    mean_other = sums / counts
    mean_other[np.arange(len(labels)), inv] = np.inf
    b = mean_other.min(axis=1)
    s = np.where(own_n > 1, (b - a) / np.maximum(np.maximum(a, b), 1e-12), 0.0)
    return float(s.mean())


def split_submodules(Cc: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """(labels 1..k, 0 = unassigned; best silhouette; leaf order) from a correlation matrix."""
    n = Cc.shape[0]
    D = np.clip(1 - Cc, 0, 2)
    np.fill_diagonal(D, 0)
    Zl = linkage(squareform(D, checks=False), method="average")
    order = leaves_list(Zl)
    best_s, best_lab = -1.0, np.ones(n, dtype=int)
    for k in range(2, min(MAX_SUBMODULES, n // MIN_SUBMODULE) + 1):
        lab = fcluster(Zl, t=k, criterion="maxclust")
        s = silhouette_precomputed(D, lab)
        if np.isfinite(s) and s > best_s:
            best_s, best_lab = s, lab
    if best_s < MIN_SILHOUETTE:
        return np.ones(n, dtype=int), best_s, order
    sizes = pd.Series(best_lab).value_counts()
    small = sizes[sizes < MIN_SUBMODULE].index
    lab = np.where(np.isin(best_lab, small), 0, best_lab)
    # renumber 1..m by size
    big = [c for c in sizes.index if c not in small]
    remap = {c: i + 1 for i, c in enumerate(big)}
    lab = np.array([remap.get(c, 0) for c in lab])
    return lab, best_s, order


def eigengene(Zz: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Mean of the members' z-scored profiles -- a robust module score per cluster."""
    return Zz[idx].mean(axis=0)


def submodule_null(ds: str, disc: str, held: str, shared: pd.Index, idx_shared: np.ndarray,
                   data: dict, rng: np.random.Generator) -> pd.DataFrame:
    """Held-out coherence of sub-modules found in random lists matched to this list.

    Any gene set splits into groups that co-vary reproducibly -- the transcriptome
    has strong, real structure (neuronal vs progenitor genes, ...). So the
    question for a list's sub-module is not "does it replicate" but "is it
    tighter, in the held-out donors, than what the same procedure carves out of
    a random list of matched genes". Random lists match the real one bin for
    bin (level x spread among genes expressed in both chemistries).
    """
    lc_d = data[(ds, disc)][0].loc[shared]
    Zd = data[(ds, disc)][2]["across_clusters"][0]
    pos_d = {g: i for i, g in enumerate(data[(ds, disc)][0].index)}
    pos_h = {g: i for i, g in enumerate(data[(ds, held)][0].index)}
    rows_d = np.array([pos_d[g] for g in shared])
    rows_h = np.array([pos_h[g] for g in shared])
    bins = C.expression_bins(lc_d)
    rand = C.matched_sets(bins, idx_shared, N_RANDOM_LISTS, rng)
    out = []
    for r in rand:
        Zr = Zd[rows_d[r]]
        lab = split_submodules(C.dot(Zr, Zr.T))[0]
        if len(set(lab) - {0}) < 2:
            continue                       # no split: nothing comparable to record
        for lv in set(lab) - {0}:
            members = rows_h[r[lab == lv]]
            rec = {"size": int(members.size)}
            for cname, (Zn, _, _, _) in data[(ds, held)][2].items():
                S = Zn[members].sum(axis=0)
                rec[cname] = float(C.coherence_from_sums(S[None, :], members.size)[0])
            out.append(rec)
    return pd.DataFrame(out)


def selection_p(null: pd.DataFrame, size: int, coh: float, cname: str) -> tuple[float, float, int]:
    """(p, null median, n) against random sub-modules of comparable size."""
    if null.empty or cname not in null:
        return np.nan, np.nan, 0
    near = null[(null["size"] >= size / 2) & (null["size"] <= size * 2)]
    if len(near) < MIN_NULL_SUBMODULES:
        near = null.iloc[(null["size"] - size).abs().argsort()[:MIN_NULL_SUBMODULES]]
    vals = near[cname].to_numpy(float)
    return float(((vals >= coh).sum() + 1) / (vals.size + 1)), float(np.median(vals)), int(vals.size)


# ---------------------------------------------------------------------------
def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    lists = C.gene_lists()
    src = C.gene_lists_label()
    if not lists:
        C.log(f"  no gene lists found in {src} -- set AIM_GENE_LISTS or copy them there")
        out.summary(TITLE, "Skipped: no gene lists.",
                    [f"Gene-list folder: `{src}`."],
                    [f"**No gene lists found in `{src}`.** Copy them there or set AIM_GENE_LISTS."],
                    [], [])
        return
    rng = np.random.default_rng(C.SEED)

    data = {}           # (ds, chem) -> (lc, annot, contexts)
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        grouping = C.COEXPR_GROUPING[ds]
        clustering = grouping.replace("cluster_", "", 1)
        out.used(f"{n}/09_pseudobulk/{grouping}__pseudobulk_counts.csv",
                 f"{n}/04_clusters/cluster_profile_{clustering}.csv")
        lc, annot = C.cluster_expression(n)
        data[(ds, chem)] = (lc, annot, C.coexpr_contexts(lc, annot))
        C.log(f"  {n}: {lc.shape[0]:,} genes x {lc.shape[1]} clusters ({grouping}, "
              f">= {C.COEXPR_MIN_CELLS} cells, donor-dominated dropped); contexts: "
              + ", ".join(data[(ds, chem)][2]))

    # ---- whole-list coherence and per-gene connectivity ---------------------
    coh_rows, gene_rows = [], []
    for (ds, chem), (lc, annot, ctx) in data.items():
        mapped = C.mapped_lists(ds)
        pos = {g: i for i, g in enumerate(lc.index)}
        for name in lists:
            idx = np.array([pos[g] for g in mapped.get(name, []) if g in pos])
            if idx.size < MIN_GENES:
                continue
            for cname, (Zn, bins, cols, _) in ctx.items():
                coh, null, conn, p_gene = coherence_test(Zn, bins, idx, rng)
                sd = null.std()
                coh_rows.append({
                    "dataset": ds, "chemistry": chem, "context": cname, "gene_list": name,
                    "n_genes": int(idx.size), "n_clusters": len(cols),
                    "coherence_mean_r": coh, "null_mean_r": float(null.mean()),
                    "null_sd": float(sd),
                    "effect_vs_null_sd": (coh - null.mean()) / sd if sd > 0 else np.nan,
                    "perm_p": float(((null >= coh).sum() + 1) / (N_RANDOM + 1))})
                for g, cval, pval in zip(lc.index[idx], conn, p_gene):
                    gene_rows.append({"dataset": ds, "chemistry": chem, "context": cname,
                                      "gene_list": name, "gene": g,
                                      "connectivity": float(cval), "perm_p": float(pval)})
    coh = pd.DataFrame(coh_rows)
    out.write(coh, "coherence_per_stratum",
              "List x stratum x context: mean pairwise Spearman r vs matched random sets")

    comb = []
    for (ds, cname, name), g in coh.groupby(["dataset", "context", "gene_list"], sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        comb.append({"dataset": ds, "context": cname, "gene_list": name,
                     "n_genes_v2": int(a.n_genes), "n_genes_v3": int(b.n_genes),
                     "coherence_v2": a.coherence_mean_r, "null_v2": a.null_mean_r,
                     "effect_v2": a.effect_vs_null_sd, "p_v2": a.perm_p,
                     "coherence_v3": b.coherence_mean_r, "null_v3": b.null_mean_r,
                     "effect_v3": b.effect_vs_null_sd, "p_v3": b.perm_p,
                     "combined_p": one_sided_combine(a.perm_p, b.perm_p)})
    comb = pd.DataFrame(comb)
    if not comb.empty:
        comb["combined_q"] = np.nan
        for _, ix in comb.groupby(["dataset", "context"]).groups.items():
            comb.loc[ix, "combined_q"] = C.bh(comb.loc[ix, "combined_p"])
        comb["tier"] = [tier(a, b, q) for a, b, q in zip(comb.p_v2, comb.p_v3, comb.combined_q)]
        comb = comb.sort_values(["dataset", "context", "combined_p"])
    out.write(comb, "coherence_combined",
              "v2 x v3 combined coherence per list and context; tier replicated / supported")

    genes = pd.DataFrame(gene_rows)
    gm = genes.pivot_table(index=["dataset", "context", "gene_list", "gene"], columns="chemistry",
                           values=["connectivity", "perm_p"]).reset_index()
    gm.columns = ["_".join(c).strip("_") for c in gm.columns]
    if {"perm_p_v2", "perm_p_v3"} <= set(gm.columns):
        gm["status"] = np.where((gm.perm_p_v2 < 0.05) & (gm.perm_p_v3 < 0.05)
                                & (gm.connectivity_v2 > 0) & (gm.connectivity_v3 > 0), "core",
                                np.where(gm[["perm_p_v2", "perm_p_v3"]].isna().any(axis=1),
                                         "one chemistry only", "peripheral"))
    out.write(gm, "gene_connectivity",
              "Per list gene: mean correlation with the other members (v2, v3) vs matched null; "
              "core = beyond null in both chemistries")

    # ---- sub-modules, cross-validated between donor sets ----------------------
    # A split chosen because its groups correlate will always "pass" a test on
    # the data it was chosen from (a random list did, in development). So each
    # chemistry discovers sub-modules and the OTHER chemistry -- disjoint donors
    # -- tests them. The consensus split is kept for display only.
    sub_rows, sub_stats, cand_rows, heat = [], [], [], []
    all_lists = {name: set(g) for name, g in lists.items()}
    other = {"v2": "v3", "v3": "v2"}
    for ds in C.DATASETS:
        mapped = C.mapped_lists(ds)
        lcs = {c: data[(ds, c)][0] for c in C.CHEMISTRIES}
        shared = lcs["v2"].index.intersection(lcs["v3"].index)
        pos = {c: {g: i for i, g in enumerate(lcs[c].index)} for c in C.CHEMISTRIES}
        for name in lists:
            members = [g for g in dict.fromkeys(mapped.get(name, [])) if g in shared]
            if len(members) < 2 * MIN_SUBMODULE:
                continue
            R = {}
            for c in C.CHEMISTRIES:
                Zn = data[(ds, c)][2]["across_clusters"][0]
                ix = np.array([pos[c][g] for g in members])
                R[c] = C.dot(Zn[ix], Zn[ix].T)
            Cc = np.tanh((np.arctanh(np.clip(R["v2"], -0.999, 0.999))
                          + np.arctanh(np.clip(R["v3"], -0.999, 0.999))) / 2)
            np.fill_diagonal(Cc, 1.0)
            lab_c, _, order = split_submodules(Cc)
            heat.append((ds, name, members, Cc, lab_c, order))
            labs = {c: split_submodules(R[c])[0] for c in C.CHEMISTRIES}
            for i, g in enumerate(members):
                sub_rows.append({"dataset": ds, "gene_list": name, "gene": g,
                                 "submodule_found_in_v2": f"{name}.v2.{labs['v2'][i]}" if labs["v2"][i] else "",
                                 "submodule_found_in_v3": f"{name}.v3.{labs['v3'][i]}" if labs["v3"][i] else "",
                                 "consensus_group_display_only": int(lab_c[i])})
            idx_shared = np.array([shared.get_loc(g) for g in members])
            for disc in C.CHEMISTRIES:
                held = other[disc]
                lab = labs[disc]
                split = len(set(lab) - {0}) >= 2
                null_sub = (submodule_null(ds, disc, held, shared, idx_shared, data, rng)
                            if split else pd.DataFrame())
                for lv in sorted(set(lab) - {0}):
                    # No split (one group, maybe a few leftovers): the unit is the whole
                    # list, which was fixed in advance, so the plain matched null applies.
                    mem = (members if not split else
                           [g for g, x in zip(members, lab) if x == lv])
                    row = {"dataset": ds, "gene_list": name, "submodule": f"{name}.{disc}.{lv}",
                           "found_in": disc, "tested_in": held, "n_genes": len(mem),
                           "no_split": not split, "genes": "|".join(mem)}
                    eig = {}
                    for chem in C.CHEMISTRIES:
                        lc, annot, ctx = data[(ds, chem)]
                        idx = np.array([pos[chem][g] for g in mem])
                        tag = "heldout" if chem == held else "insample"
                        for cname, (Zn, bins, cols, Zz) in ctx.items():
                            c_obs, null, _, _ = coherence_test(Zn, bins, idx, rng)
                            row[f"coherence_{cname}_{tag}"] = c_obs
                            row[f"null_{cname}_{tag}"] = float(null.mean())
                            row[f"p_{cname}_{tag}"] = float(((null >= c_obs).sum() + 1) / (N_RANDOM + 1))
                        Zn, bins, cols, Zz = ctx["across_clusters"]
                        e = eigengene(Zz, idx)
                        eig[chem] = (e, Zn, lc.index)
                        a = annot.loc[cols]
                        by_cls = pd.Series(e, index=cols).groupby(a["dominant_cell_class"].to_numpy()).mean()
                        by_cls = by_cls[a["dominant_cell_class"].value_counts().reindex(by_cls.index) >= 3]
                        top = pd.Series(e, index=cols).sort_values(ascending=False).head(3)
                        row[f"peak_class_{chem}"] = by_cls.idxmax() if len(by_cls) else ""
                        row[f"top_clusters_{chem}"] = "; ".join(
                            f"{c} ({a.loc[c, 'dominant_cell_class']}, {a.loc[c, 'dominant_region']}, "
                            f"{a.loc[c, 'age_pcw_median']:g} pcw)" for c in top.index)
                        row[f"rho_with_cluster_age_{chem}"] = float(
                            stats.spearmanr(e, a["age_pcw_median"].to_numpy(float)).statistic)
                    # significance: whole list -> plain matched null (the list is fixed in
                    # advance); a split-off group -> selection-aware null from random lists
                    for cname in CONTEXTS:
                        if f"coherence_{cname}_heldout" not in row:
                            continue
                        if row["no_split"]:
                            row[f"p_select_{cname}"] = row[f"p_{cname}_heldout"]
                            row[f"random_submodule_median_{cname}"] = row[f"null_{cname}_heldout"]
                        else:
                            pv, med, nn = selection_p(null_sub, len(mem),
                                                      row[f"coherence_{cname}_heldout"], cname)
                            row[f"p_select_{cname}"] = pv
                            row[f"random_submodule_median_{cname}"] = med
                            row[f"n_random_submodules_{cname}"] = nn
                    sub_stats.append(row)
                    rs = {}
                    for chem, (e, Zn, idxs) in eig.items():
                        ez = C.unit_rank_rows(e[None, :])[0]
                        rs[chem] = pd.Series(C.dot(Zn, ez), index=idxs)
                    both = rs["v2"].index.intersection(rs["v3"].index)
                    cand = pd.DataFrame({"r_v2": rs["v2"].loc[both], "r_v3": rs["v3"].loc[both]})
                    cand = cand[~cand.index.isin(mem)]
                    cand["r_min"] = cand.min(axis=1)
                    cand = cand[cand.r_min >= CAND_R].sort_values("r_min", ascending=False).head(N_CAND)
                    for g, r in cand.iterrows():
                        cand_rows.append({"dataset": ds, "submodule": row["submodule"], "gene": g,
                                          "r_v2": r.r_v2, "r_v3": r.r_v3, "r_min": r.r_min,
                                          "in_other_lists": "|".join(k for k, v in all_lists.items()
                                                                     if g in v and k != name)})
    subs = pd.DataFrame(sub_rows)
    sstats = pd.DataFrame(sub_stats)
    if not sstats.empty:
        for cname in CONTEXTS:
            col = f"p_select_{cname}"
            if col not in sstats:
                continue
            sstats[f"q_select_{cname}"] = np.nan
            for _, ix in sstats.groupby("dataset").groups.items():
                ok = sstats.loc[ix, col].notna()
                sstats.loc[ok[ok].index, f"q_select_{cname}"] = C.bh(sstats.loc[ok[ok].index, col])
        sstats["validated"] = sstats["q_select_across_clusters"] < 0.05
        # Pair each sub-module with its best match discovered in the other chemistry.
        sstats["best_match"], sstats["jaccard"], sstats["found_both_ways"] = "", np.nan, False
        sstats["core_genes"] = ""
        for i, r in sstats.iterrows():
            cand_m = sstats[(sstats.dataset == r.dataset) & (sstats.gene_list == r.gene_list)
                            & (sstats.found_in != r.found_in)]
            if cand_m.empty:
                continue
            A = set(r.genes.split("|"))
            jac = cand_m.genes.map(lambda x: len(A & set(x.split("|"))) / len(A | set(x.split("|"))))
            j = jac.idxmax()
            B = set(cand_m.loc[j, "genes"].split("|"))
            sstats.loc[i, ["best_match", "jaccard"]] = [cand_m.loc[j, "submodule"], jac.max()]
            sstats.loc[i, "found_both_ways"] = bool(r.validated and cand_m.loc[j, "validated"]
                                                    and jac.max() >= 0.3)
            sstats.loc[i, "core_genes"] = "|".join(g for g in r.genes.split("|") if g in B)
    cands = pd.DataFrame(cand_rows, columns=["dataset", "submodule", "gene", "r_v2", "r_v3",
                                             "r_min", "in_other_lists"])
    if not cands.empty and not sstats.empty:
        # sub-module names repeat across datasets: filter on (dataset, sub-module)
        ok = set(zip(sstats.loc[sstats.validated, "dataset"], sstats.loc[sstats.validated, "submodule"]))
        cands = cands[[k in ok for k in zip(cands.dataset, cands.submodule)]]
    out.write(subs, "submodule_membership",
              "List gene -> sub-module found in v2 and in v3 (and the display-only consensus group)")
    out.write(sstats, "submodules",
              "Per sub-module: found in one chemistry, coherence re-tested in the other (held "
              "out) against sub-modules carved the same way from matched random lists; best "
              "match among the other chemistry's sub-modules; peak class, top clusters")
    out.write(cands, "candidate_members",
              f"Non-members correlating >= {CAND_R} with a validated sub-module's eigengene in both chemistries")

    figures(out, comb, heat)

    # ---- findings -------------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        for cname in CONTEXTS:
            g = comb[(comb.dataset == ds) & (comb.context == cname)] if not comb.empty else comb
            if g.empty:
                continue
            rep = g[g.tier == "replicated"]
            what = ("co-expressed across clusters" if cname == "across_clusters"
                    else "co-vary within cell classes")
            f.append(f"**{ds}: lists whose genes {what}, replicated in both donor sets** -- "
                     + ("; ".join(f"{r.gene_list} (mean r {r.coherence_v2:.2f}/{r.coherence_v3:.2f} "
                                  f"vs {r.null_v2:.2f}/{r.null_v3:.2f} for matched random genes)"
                                  for r in rep.itertuples()) if len(rep) else "none")
                     + ". Not coherent: "
                     + (", ".join(g[g.tier == ""].gene_list) or "none") + ".")
    if not sstats.empty:
        for ds in C.DATASETS:
            s_ds = sstats[(sstats.dataset == ds) & sstats.validated & ~sstats.no_split]
            # one line per pair found both ways (report the v2-discovered side), then singles
            both = s_ds[s_ds.found_both_ways & (s_ds.found_in == "v2")]
            single = s_ds[~s_ds.found_both_ways]
            if both.empty and single.empty:
                n_tested = int(((sstats.dataset == ds) & ~sstats.no_split).sum())
                f.append(f"**{ds}: no sub-module beats what the same procedure finds in matched "
                         f"random lists** ({n_tested} tested); where lists are coherent, they "
                         "behave as one group.")
                continue
            parts = []
            for r in both.head(MAX_LISTED).itertuples():
                core = r.core_genes.split("|")
                parts.append(f"{r.gene_list}: {len(core)} core genes found in both donor sets "
                             f"({', '.join(core[:6])}{', ...' if len(core) > 6 else ''}; peaks in "
                             f"{r.peak_class_v2}"
                             + ("" if r.peak_class_v2 == r.peak_class_v3 else f" / {r.peak_class_v3}")
                             + (", also co-varies within classes"
                                if getattr(r, "q_select_within_class", 1) < 0.05 else "") + ")")
            f.append(f"**{ds}: sub-modules that hold up across donor sets** (found in one chemistry, "
                     "coherent in the other, and found again there) -- "
                     + ("; ".join(parts) if parts else "none")
                     + (". Validated in one direction only: "
                        + "; ".join(f"{r.submodule} ({r.n_genes} genes, peaks in {r.peak_class_v2})"
                                    for r in single.head(MAX_LISTED).itertuples())
                        if len(single) else "") + ".")
    if not cands.empty:
        # A pair found both ways shares its candidates; list it once (its v2 side).
        dup = set(zip(sstats.loc[sstats.found_both_ways & (sstats.found_in == "v3"), "dataset"],
                      sstats.loc[sstats.found_both_ways & (sstats.found_in == "v3"), "submodule"]))
        shown = cands[[k not in dup for k in zip(cands.dataset, cands.submodule)]]
        top = (shown.sort_values("r_min", ascending=False)
               .groupby(["dataset", "submodule"]).head(5)
               .groupby(["dataset", "submodule"])["gene"].apply(", ".join))
        f.append("**Candidate new members** (correlate >= "
                 f"{CAND_R} with a validated sub-module in both chemistries; hypotheses): "
                 + "; ".join(f"{sm} ({ds}): {genes}" for (ds, sm), genes in top.head(MAX_LISTED).items())
                 + ("" if len(top) <= MAX_LISTED else f" (+{len(top) - MAX_LISTED} more)")
                 + ". Full list: candidate_members.csv.")

    out.summary(
        TITLE,
        "Across the fine clusters of each file, do a list's genes rise and fall together "
        "beyond what genes of the same expression level and spread do; which members carry "
        "that; does the list split into tighter sub-modules; and which outside genes track them?",
        [f"Clusters: {C.COEXPR_GROUPING['cortex']} (cortex) and {C.COEXPR_GROUPING['human_dev']} "
         f"(human_dev), >= {C.COEXPR_MIN_CELLS} cells, donor-dominated clusters dropped; log2 "
         f"TMM-CPM; genes >= {C.COEXPR_MIN_CPM:g} CPM in >= {C.COEXPR_MIN_CLUSTERS} clusters; "
         "Spearman correlation.",
         "Contexts: across all clusters; and within cell classes (each gene's mean per "
         f"dominant class removed; classes with >= {C.COEXPR_MIN_CLUSTERS_PER_CLASS} clusters).",
         f"Coherence = mean pairwise r; null = {N_RANDOM:,} random sets matched member by "
         "member on mean level x spread (5 x 5 quantile bins), no repeats; one-sided p; v2 and v3 "
         "combined by Stouffer, BH per dataset x context; replicated = both nominal + q < 0.05.",
         "Sub-modules: each chemistry splits a list on its own correlations (average linkage, "
         f"best silhouette, 2-{MAX_SUBMODULES} groups, none below {MIN_SILHOUETTE}; groups < "
         f"{MIN_SUBMODULE} genes unassigned); each sub-module is re-tested in the other "
         f"chemistry, against sub-modules carved the same way from {N_RANDOM_LISTS} random lists "
         "matched to the real one (held-out coherence, size within 2-fold; BH per dataset). "
         "In-sample coherence is reported but never used for significance. The heatmaps show "
         "the Fisher-z average of both chemistries, for display only.",
         "Eigengene = mean z-scored member profile per cluster; peak class = class with the "
         "highest mean eigengene (classes with >= 3 clusters)."],
        f,
        ["Cluster pseudobulks are not independent samples (clusters nest in classes, donors, "
         "regions); the matched null asks whether the list beats comparable genes in the same "
         "matrix, not whether correlations are 'significant' in absolute terms.",
         "Across-cluster coherence is largely cell-identity co-expression; the within-class "
         "context is the stricter test of co-regulation.",
         "Null sets ignore gene-gene correlation structure beyond level and spread; lists of "
         "paralogs or physically clustered genes can look coherent for that reason.",
         "cortex and human_dev share donors; agreement between them is not replication."],
        ["Test sub-modules for shared regulators (TF motif or ChIP target enrichment).",
         "Repeat within single regions of human_dev once a region x cluster export exists."])


def figures(out: C.Output, comb: pd.DataFrame, heat: list) -> None:
    plt = C.plt_or_none()
    if plt is None:
        return
    if not comb.empty:
        for ds in C.DATASETS:
            g = comb[comb.dataset == ds]
            if g.empty:
                continue
            names = sorted(g.gene_list.unique())
            fig, axes = plt.subplots(1, len(CONTEXTS), figsize=(5 * len(CONTEXTS), 0.4 * len(names) + 1.8),
                                     squeeze=False)
            for ax, cname in zip(axes.flat, CONTEXTS):
                h = g[g.context == cname].set_index("gene_list").reindex(names)
                y = np.arange(len(names))
                ax.barh(y - 0.2, h.effect_v2, height=0.4, label="v2")
                ax.barh(y + 0.2, h.effect_v3, height=0.4, label="v3")
                ax.axvline(0, color="grey", lw=0.8)
                ax.set_yticks(y, names, fontsize=7)
                ax.set_xlabel("coherence vs matched random genes (null SDs)", fontsize=8)
                ax.set_title(f"{ds}: {cname.replace('_', ' ')}", fontsize=9)
            axes.flat[0].legend(fontsize=7)
            fig.tight_layout()
            out.figure(fig, f"coherence_{ds}", f"List coherence per chemistry and context, {ds}")
            plt.close(fig)
    for ds, name, members, Cc, lab, order in heat:
        if len(members) > HEATMAP_MAX_GENES:
            keep = np.argsort(-np.abs(Cc).mean(axis=1))[:HEATMAP_MAX_GENES]
            order = [i for i in order if i in set(keep)]
        M = Cc[np.ix_(order, order)]
        fig, ax = plt.subplots(figsize=(7, 6.4))
        im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
        ticks = [members[i] for i in order]
        if len(ticks) <= 80:
            ax.set_xticks(range(len(ticks)), ticks, rotation=90, fontsize=5)
            ax.set_yticks(range(len(ticks)), ticks, fontsize=5)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        labs = np.array(lab)[order]
        for b in np.nonzero(labs[1:] != labs[:-1])[0]:
            ax.axhline(b + 0.5, color="k", lw=0.6)
            ax.axvline(b + 0.5, color="k", lw=0.6)
        ax.set_title(f"{ds}: {name} -- co-expression, mean of v2/v3 (boxes: display-only grouping)",
                     fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.7, label="Spearman r (Fisher-z mean of v2, v3)")
        fig.tight_layout()
        out.figure(fig, f"coexpression_{ds}_{name}", f"Consensus co-expression of {name}, {ds}")
        plt.close(fig)


if __name__ == "__main__":
    main()
