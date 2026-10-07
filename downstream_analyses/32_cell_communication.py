#!/usr/bin/env python3
"""32 - Ligand-receptor signalling between cell types over development.

Question: which cell types could signal to which (a ligand expressed by the
sender, its receptor by the receiver), how does that potential change with
age in both donor sets, and where do NDD genes act as ligands or receptors?

Method
  Pairs: OmniPath interactions whose source is a ligand (secreted or on the
  membrane) and target a receptor (plasma membrane), each so annotated by >= 4
  intercell resources (fewer admits e.g. PCNA as a 'ligand'); symbols mapped to each dataset's.
  Expression: pseudobulk counts per cell class x age, TMM-normalised jointly
  over all classes and ages of a stratum (so classes are comparable), as CPM.
  Score of pair (L, R) from sender s to receiver r at an age = sqrt(CPM of L
  in s x CPM of R in r); the pair is 'expressed' there when both reach 10 CPM.
  Per stratum and (pair, sender, receiver) expressed at >= half of the >= 5
  age points both classes share: Spearman of log2(score + 1) with age, exact
  permutation (shared null per sender-receiver age set); v2 x v3 signed
  Stouffer, BH per dataset, tiered.
  NDD genes: list and seed-panel genes that are ligands or receptors, the pairs
  they take part in and their trends.

Inputs: <annotations>/lr/{omnipath_interactions.tsv, intercell_ligand_receptor.tsv};
csv_exports/<ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,group_summary}.csv; gene lists
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "32_cell_communication"
TITLE = "Ligand-receptor signalling between cell types over development"
MIN_CPM = 10.0
MIN_AGES = 5
MAX_LISTED = 10


def stratum(n: str, ds: str, pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cnt = C.group_matrix(n, "cell_class_x_age", "pseudobulk_counts")
    lc = C.tmm_log_cpm(cnt)
    cpm = 2 ** lc - 1
    lab = [C.split_class_age(c) for c in cpm.columns]
    keep = [i for i, (_, a) in enumerate(lab) if not C.excluded(ds, age=a)]
    cpm = cpm.iloc[:, keep]
    lab = [lab[i] for i in keep]
    classes = sorted({c for c, _ in lab})
    ages = sorted({a for _, a in lab})
    col = {(c, a): j for j, (c, a) in enumerate(lab)}
    p = pairs[pairs.ligand.isin(cpm.index) & pairs.receptor.isin(cpm.index)].reset_index(drop=True)
    L = cpm.loc[p.ligand].to_numpy(float)
    R = cpm.loc[p.receptor].to_numpy(float)
    rows, overview = [], []
    for s in classes:
        for r in classes:
            shared = [a for a in ages if (s, a) in col and (r, a) in col]
            if len(shared) < MIN_AGES:
                continue
            Ls = L[:, [col[(s, a)] for a in shared]]
            Rr = R[:, [col[(r, a)] for a in shared]]
            expr = (Ls >= MIN_CPM) & (Rr >= MIN_CPM)
            overview.append({"sender": s, "receiver": r, "n_ages": len(shared),
                             "mean_expressed_pairs": float(expr.sum(axis=0).mean())})
            ok = expr.sum(axis=1) >= len(shared) / 2
            if not ok.any():
                continue
            score = np.log2(np.sqrt(Ls[ok] * Rr[ok]) + 1)
            x = np.array(shared)
            rho = C.spearman_rows(score, x)
            pp, _ = C.spearman_perm_p(rho, x)
            sub = p[ok]
            rows.append(pd.DataFrame({"ligand": sub.ligand.to_numpy(), "receptor": sub.receptor.to_numpy(),
                                      "sender": s, "receiver": r, "n_ages": len(shared), "rho_vs_age": rho,
                                      "perm_p": pp, "mean_log2_score": score.mean(axis=1)}))
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(), pd.DataFrame(overview))


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    lr = C.ligand_receptor_pairs()
    if lr.empty:
        out.summary(TITLE, "Skipped: ligand-receptor annotations missing.",
                    [f"Needs lr/omnipath_interactions.tsv and lr/intercell_ligand_receptor.tsv in {C.annotations_dir()}."],
                    ["**Not run**: fetch the annotation files first (running_scripts/fetch_annotations.sh), then rerun "
                     "step 3."], [], [])
        return
    out.used("annotations/lr/omnipath_interactions.tsv", "annotations/lr/intercell_ligand_receptor.tsv")
    per, ov, ndd_rows = [], [], []
    for ds in C.DATASETS:
        mp = C.to_dataset_symbols(sorted(set(lr.ligand) | set(lr.receptor)), ds)
        pairs = lr.assign(ligand=lr.ligand.map(mp), receptor=lr.receptor.map(mp)).dropna()
        pairs = pairs[(pairs.ligand != "") & (pairs.receptor != "")].drop_duplicates()
        for chem in C.CHEMISTRIES:
            n = C.ns(ds, chem)
            out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv")
            r, o = stratum(n, ds, pairs)
            if len(r):
                per.append(r.assign(dataset=ds, chemistry=chem))
            if len(o):
                ov.append(o.assign(dataset=ds, chemistry=chem))
        sets = C.analysis_gene_sets(ds, groups=("ndd",), modules=False, out=out)
        for name, genes in sets.items():
            g = set(genes)
            for role, col in (("ligand", "ligand"), ("receptor", "receptor")):
                hit = pairs[pairs[col].isin(g)]
                for gene, h in hit.groupby(col):
                    ndd_rows.append({"dataset": ds, "gene_set": name, "gene": gene, "role": role,
                                     "partners": "|".join(sorted(set(h.receptor if role == "ligand" else h.ligand))[:20])})
    if not per:
        out.summary(TITLE, "No sender-receiver pair has enough shared age points.", [], ["**Nothing testable.**"], [], [])
        return
    per = pd.concat(per, ignore_index=True)
    ov = pd.concat(ov, ignore_index=True)
    ndd = pd.DataFrame(ndd_rows)
    comb = C.combine_chemistries(per, ["dataset", "ligand", "receptor", "sender", "receiver"], effect="rho_vs_age",
                                 weight="n_ages", labels=("rises with age", "falls with age"),
                                 carry=("rho_vs_age", "mean_log2_score"))
    ndd_genes = set(ndd.gene) if len(ndd) else set()
    if not comb.empty:
        comb["ndd_gene"] = comb.ligand.isin(ndd_genes) | comb.receptor.isin(ndd_genes)
    out.write(ov, "sender_receiver_overview", f"Per stratum x sender x receiver: mean number of pairs with both partners "
              f">= {MIN_CPM:g} CPM per age point")
    out.write(per, "interaction_age_trends_per_stratum", "Per stratum x pair x sender x receiver: Spearman of log2 score with age")
    out.write(comb, "interaction_age_trends_combined", "v2 x v3 combined; BH per dataset; tier; NDD gene flag")
    out.write(ndd, "ndd_ligands_receptors", "NDD list / seed-panel genes that are ligands or receptors, with partners")

    f = []
    o = ov.groupby(["dataset", "sender", "receiver"]).mean_expressed_pairs.mean().reset_index()
    for ds, g in o.groupby("dataset"):
        top = g.sort_values("mean_expressed_pairs", ascending=False).head(6)
        f.append(f"**{ds}: busiest sender -> receiver routes** (pairs expressed per age, mean of chemistries): "
                 + ", ".join(f"{r.sender} -> {r.receiver} {r.mean_expressed_pairs:.0f}" for r in top.itertuples()) + ".")
    if not comb.empty:
        for ds in C.DATASETS:
            g = comb[(comb.dataset == ds) & (comb.tier == "replicated")]
            tested = int((comb.dataset == ds).sum())
            routes = g.groupby(["sender", "receiver", "direction"]).size().unstack(fill_value=0) if len(g) else pd.DataFrame()
            f.append(f"**{ds}: replicated changes in signalling potential** ({len(g)} of {tested} pair x route tests): "
                     + ("; ".join(f"{s} -> {r}: {int(row.get('rises with age', 0))} rise / {int(row.get('falls with age', 0))} fall"
                                  for (s, r), row in routes.sort_values(list(routes.columns)[0], ascending=False).head(8).iterrows())
                        if len(routes) else "none") + ".")
            nd = g[g.ndd_gene].sort_values("combined_p")
            if len(nd):
                f.append(f"**{ds}: replicated changes involving NDD genes** (ligand -> receptor, sender -> receiver; "
                         "rho v2/v3): " + "; ".join(
                             f"{r.ligand} -> {r.receptor} ({r.sender} -> {r.receiver}) {r.direction} "
                             f"({r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f})" for r in nd.head(MAX_LISTED).itertuples())
                         + ("" if len(nd) <= MAX_LISTED else f" (+{len(nd) - MAX_LISTED} more)") + ".")
    if len(ndd):
        lists = ndd[ndd.gene_set.str.startswith("list:")]
        if len(lists):
            f.append("**NDD list genes that are ligands or receptors**: " + "; ".join(
                f"{ds} {gs}: " + ", ".join(sorted(set(h.gene))[:12]) + (" ..." if h.gene.nunique() > 12 else "")
                for (ds, gs), h in lists.groupby(["dataset", "gene_set"])) + ".")
    out.summary(
        TITLE,
        "Which cell types could signal to which through ligand-receptor pairs, how does that potential change "
        "with age in both donor sets, and where do NDD genes act as ligands or receptors?",
        ["Pairs: OmniPath interactions, ligand (secreted / membrane) -> receptor (plasma membrane), each role "
         "annotated by >= 4 intercell resources; symbols mapped to each dataset.",
         f"Score = sqrt(ligand CPM in sender x receptor CPM in receiver), CPM from TMM over all class x age "
         f"pseudobulks of a stratum; expressed = both >= {MIN_CPM:g} CPM.",
         f"Age trend per pair x route expressed at >= half of >= {MIN_AGES} shared age points: Spearman, exact "
         "permutation; v2 x v3 signed Stouffer, BH per dataset, tiered."],
        f,
        ["Expression of a ligand and its receptor is potential, not signalling: protein levels, processing, "
         "localisation and spatial proximity are not measured.",
         "Complex receptors (multiple subunits) are reduced to pairwise interactions.",
         "Within a class, a trend can come from sub-type mix (as everywhere in 03)."],
        ["Spatial data to check which sender-receiver pairs are neighbours."])


if __name__ == "__main__":
    main()
