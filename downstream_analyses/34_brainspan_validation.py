#!/usr/bin/env python3
"""34 - Do the cortex age trends replicate in an independent atlas (BrainSpan)?

Question: cortex donors are a subset of human_dev's, so the two files cannot
replicate each other. BrainSpan's developmental transcriptome (bulk RNA-seq of
dissected brain regions, other donors, other lab) can. Do genes and gene lists
that change with age in cortex change the same way in BrainSpan neocortex over
the overlapping ages?

Method
  BrainSpan (RNA-seq, Gencode v10 genes, RPKM): neocortical samples (all
  cortical areas) of donors aged 8-16 pcw (our cortex spans 5.5-14 pcw), log2
  (RPKM + 1) averaged per donor; genes >= 1 RPKM on average; Spearman with donor
  age, permutation p. Genes matched to cortex symbols via Ensembl id.
  Ours: whole-cortex pseudobulk per age point (all cells -- comparable to bulk
  tissue), log2 TMM-CPM, genes >= 5 CPM, Spearman with age per chemistry, v2 x v3
  signed Stouffer.
  Agreement: Spearman of our combined Z with BrainSpan's rho over shared genes;
  for our replicated genes, the share with the same sign in BrainSpan (50% by
  chance), also among those nominal in BrainSpan. Lists and seed panels: mean
  BrainSpan rho of members vs random genes matched on BrainSpan expression
  decile, beside our whole-cortex list trend.

Inputs: <annotations>/brainspan/{expression_matrix,rows_metadata,columns_metadata}.csv;
csv_exports/cortex__<chem>/09_pseudobulk/age__{pseudobulk_counts,group_summary}.csv; gene lists
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "34_brainspan_validation"
TITLE = "External validation: cortex age trends vs BrainSpan neocortex"
DATASET = "cortex"
AGE_MIN, AGE_MAX = 8.0, 16.0
NEOCORTEX = {"Ocx", "M1C-S1C", "MFC", "DFC", "VFC", "OFC", "PCx", "TCx", "IPC", "A1C", "ITC", "STC", "V1C",
             "S1C", "M1C"}
MIN_RPKM = 1.0
MIN_CPM = 5.0
N_RANDOM = 2000
MAX_LISTED = 10


def brainspan() -> tuple[pd.DataFrame, pd.DataFrame] | None:
    base = C.annotations_dir() / "brainspan"
    if not (base / "expression_matrix.csv").exists():
        return None
    cols = pd.read_csv(base / "columns_metadata.csv")
    rows = pd.read_csv(base / "rows_metadata.csv")
    cols["pcw"] = pd.to_numeric(cols.age.str.extract(r"^([\d.]+) pcw")[0], errors="coerce")
    keep = cols.pcw.between(AGE_MIN, AGE_MAX) & cols.structure_acronym.isin(NEOCORTEX)
    use = cols.index[keep].to_numpy()
    X = pd.read_csv(base / "expression_matrix.csv", header=None, index_col=0,
                    usecols=[0] + list(use + 1)).to_numpy(float)
    L = np.log2(X + 1)
    donors = cols.loc[keep, "donor_id"].to_numpy()
    ud = list(dict.fromkeys(donors))
    M = np.stack([L[:, donors == d].mean(axis=1) for d in ud], axis=1)
    ages = cols.loc[keep].groupby("donor_id").pcw.first().reindex(ud).to_numpy()
    expr = pd.DataFrame(M, index=rows.ensembl_gene_id.astype(str).to_numpy(), columns=[str(d) for d in ud])
    meta = pd.DataFrame({"donor": [str(d) for d in ud], "age_pcw": ages,
                         "n_samples": [int((donors == d).sum()) for d in ud]})
    rpkm_mean = (2 ** expr - 1).mean(axis=1)
    expr = expr[rpkm_mean >= MIN_RPKM]
    expr = expr[~expr.index.duplicated()]
    return expr, meta


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    bs = brainspan()
    if bs is None:
        out.summary(TITLE, "Skipped: BrainSpan not fetched.",
                    [f"Needs brainspan/expression_matrix.csv in {C.annotations_dir()}."],
                    ["**Not run**: fetch the annotation files first (running_scripts/fetch_annotations.sh), then rerun "
                     "step 3."], [], [])
        return
    expr, meta = bs
    out.used("annotations/brainspan/expression_matrix.csv", "annotations/brainspan/rows_metadata.csv",
             "annotations/brainspan/columns_metadata.csv")
    C.log(f"  BrainSpan: {len(meta)} donors {meta.age_pcw.min():g}-{meta.age_pcw.max():g} pcw, {len(expr):,} genes")
    mp = C.to_dataset_symbols(list(expr.index), DATASET)
    expr = expr[expr.index.isin(mp)]
    expr.index = [mp[g] for g in expr.index]
    expr = expr[~expr.index.duplicated()]
    ages = meta.age_pcw.to_numpy(float)
    bt = C.gene_age_trends(expr, ages).set_index("gene")
    bt["bs_mean_log2_rpkm"] = expr.mean(axis=1)

    per = []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        out.used(f"{n}/09_pseudobulk/age__pseudobulk_counts.csv")
        cnt = C.group_matrix(n, "age", "pseudobulk_counts")
        a = np.array([float(c) for c in cnt.columns])
        lc = C.tmm_log_cpm(cnt)
        lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
        t = C.gene_age_trends(lc, a).assign(chemistry=chem, dataset=DATASET, n_ages=len(a))
        per.append(t)
    per = pd.concat(per, ignore_index=True)
    ours = C.combine_chemistries(per, ["dataset", "gene"], effect="rho", weight="n_ages",
                                 labels=("rises with age", "falls with age"), carry=("rho",)).set_index("gene")
    shared = ours.index.intersection(bt.index)
    j = ours.loc[shared].join(bt[["rho", "perm_p", "bs_mean_log2_rpkm"]].rename(
        columns={"rho": "brainspan_rho", "perm_p": "brainspan_p"}))
    rep = j[j.tier == "replicated"]
    same = np.sign(rep.stouffer_z) == np.sign(rep.brainspan_rho)
    nom = rep[rep.brainspan_p < 0.05]
    same_nom = np.sign(nom.stouffer_z) == np.sign(nom.brainspan_rho)
    agree = pd.DataFrame([{"n_shared_genes": len(j),
                           "spearman_ours_z_vs_brainspan_rho": float(j.stouffer_z.corr(j.brainspan_rho, method="spearman")),
                           "n_replicated_ours": len(rep), "same_sign_in_brainspan": float(same.mean()) if len(rep) else np.nan,
                           "n_replicated_and_brainspan_nominal": len(nom),
                           "same_sign_among_nominal": float(same_nom.mean()) if len(nom) else np.nan,
                           "brainspan_donors": len(meta), "brainspan_ages": f"{ages.min():g}-{ages.max():g}"}])
    sets = C.analysis_gene_sets(DATASET, groups=("ndd", "cell_cycle"), modules=True, out=out)
    bst = C.set_shift_test(bt.rho, bt.bs_mean_log2_rpkm, sets, rng, N_RANDOM)
    ost = C.set_shift_test((ours.rho_v2 + ours.rho_v3) / 2, per.groupby("gene").mean_log2cpm.mean(),
                           sets, rng, N_RANDOM)
    lists = bst.merge(ost[["gene_set", "mean", "effect_vs_null_sd", "perm_p"]], on="gene_set",
                      suffixes=("_brainspan", "_ours"), how="left") if len(bst) else pd.DataFrame()
    if len(lists):
        lists["q_brainspan"] = C.bh(lists.perm_p_brainspan)
        lists["same_direction"] = np.sign(lists.mean_brainspan) == np.sign(lists.mean_ours)
    out.write(meta, "brainspan_donors", "BrainSpan donors used: age and neocortical samples averaged")
    out.write(j.reset_index(), "gene_trends_ours_vs_brainspan", "Per gene: our whole-cortex combined trend and "
              "BrainSpan neocortex rho with age")
    out.write(agree, "agreement", "Gene-level agreement between our cortex trends and BrainSpan")
    out.write(lists, "set_trends_ours_vs_brainspan", "Per gene set: mean rho with age in BrainSpan (vs matched null) "
              "and in our whole-cortex pseudobulk")

    a = agree.iloc[0]
    f = [f"**BrainSpan**: {int(a.brainspan_donors)} donors, {a.brainspan_ages} pcw, neocortical samples averaged per "
         f"donor; {int(a.n_shared_genes):,} genes shared with cortex.",
         f"**Gene-level agreement**: Spearman of our combined Z with BrainSpan rho {a.spearman_ours_z_vs_brainspan_rho:+.2f}; "
         f"of {int(a.n_replicated_ours)} genes with a replicated whole-cortex age trend, "
         f"{a.same_sign_in_brainspan:.0%} change the same way in BrainSpan (50% by chance), "
         f"{a.same_sign_among_nominal:.0%} of the {int(a.n_replicated_and_brainspan_nominal)} also nominal there."]
    if len(lists):
        sig = lists[(lists.q_brainspan < 0.05)]
        f.append("**Gene sets moving with age in BrainSpan** (mean rho BrainSpan vs ours; same direction?): " + ("; ".join(
            f"{r.gene_set} {r.mean_brainspan:+.2f} vs {r.mean_ours:+.2f} ({'yes' if r.same_direction else 'NO'})"
            for r in sig.sort_values("perm_p_brainspan").head(MAX_LISTED * 2).itertuples()) if len(sig) else "none at q < 0.05") + ".")
    out.summary(
        TITLE,
        "Do genes and gene lists that change with age in cortex change the same way in an independent atlas "
        "(BrainSpan neocortex) over the overlapping ages?",
        [f"BrainSpan RNA-seq (Gencode v10), neocortical samples of donors {AGE_MIN:g}-{AGE_MAX:g} pcw, log2(RPKM + 1) "
         f"averaged per donor, genes >= {MIN_RPKM:g} RPKM; Spearman with age, permutation p; matched via Ensembl id.",
         f"Ours: whole-cortex pseudobulk per age point, log2 TMM-CPM, genes >= {MIN_CPM:g} CPM; Spearman per "
         "chemistry, v2 x v3 signed Stouffer.",
         f"Lists: mean BrainSpan rho vs {N_RANDOM:,} random genes matched on BrainSpan expression decile; BH."],
        f,
        ["Bulk tissue mixes cell types, and so does our whole-cortex pseudobulk: agreement here includes shifts "
         "in cell-type composition with age, which is the point of a like-for-like comparison but not a "
         "within-cell-type validation.",
         "BrainSpan has few prenatal donors in this window (about one per age point) and starts at 8 pcw; our "
         "cortex reaches back to 5.5 pcw.",
         "BrainSpan dissections are cortical areas, ours are EMX1-lineage cells; interneurons and non-neural "
         "cells are in BrainSpan only."],
        ["Deconvolve BrainSpan with our cell-class profiles to compare within-class trends."])


if __name__ == "__main__":
    main()
