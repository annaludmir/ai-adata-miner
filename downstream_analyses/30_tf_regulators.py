#!/usr/bin/env python3
"""30 - Which transcription factors are in, and may drive, the modules and NDD lists?

Question: are transcription factors (TFs) over-represented among a module's or
list's genes; whose known targets are over-represented in it (a candidate
upstream regulator); and is that TF's expression actually tied to its targets'
across clusters in these data, in both donor sets?

Method
  TFs: Lambert et al. 2018 (1,639 human TFs). Targets: CollecTRI (via
  OmniPath), signed where known. Gene symbols mapped to each dataset's.
  A. TF content: TFs among a query group's expressed genes vs the expressed
     genome (hypergeometric).
  B. Regulon enrichment: for each TF with >= 10 expressed targets,
     hypergeometric enrichment of its targets in the query group; BH within the
     group. Query groups: robust 08 modules, user lists, seed NDD panels.
  C. Regulon activity: per stratum, across the 07/08 cluster pseudobulks, the
     Spearman correlation between the TF's expression and its targets' signed
     mean Z (activated targets +, repressed -), against random target sets
     drawn from the same mean-level x spread bins; v2 x v3 combined, tiered. A
     regulator from B whose regulon is active in C is a candidate driver.
     Each TF's own within-class age trend (03) is attached.

Inputs: <annotations>/tf/{TF_names_v_1.01.txt, collectri.tsv}; results/03, 08;
gene lists; csv_exports/<ds>__<chem>/09_pseudobulk/<finest clustering> and cell_class pseudobulks
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

import _common as C

SLUG = "30_tf_regulators"
TITLE = "Transcription factors: content, target enrichment and regulon activity in modules and NDD lists"
MIN_TARGETS = 10
MIN_QUERY = 10
N_RANDOM = 200
MIN_EXPR_CPM = 1.0
MAX_LISTED = 8


def query_groups(out: C.Output, ds: str) -> dict[str, list[str]]:
    sets = C.analysis_gene_sets(ds, groups=("ndd",), modules=True, out=out)
    return sets


def regulons(ds: str, universe: set[str]) -> dict[str, tuple[list[str], list[int]]]:
    net = C.collectri()
    mp = C.to_dataset_symbols(sorted(set(net.tf) | set(net.target)), ds)
    net = net.assign(tf=net.tf.map(mp), target=net.target.map(mp)).dropna()
    net = net[(net.tf != "") & (net.target != "") & net.target.isin(universe) & (net.tf != net.target)]
    out = {}
    for tf, g in net.groupby("tf"):
        if len(g) >= MIN_TARGETS:
            out[tf] = (list(g.target), [int(s) if s != 0 else 1 for s in g.sign])
    return out


def activity(out: C.Output, ds: str, regs: dict, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for chem in C.CHEMISTRIES:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/{C.COEXPR_GROUPING[ds]}__pseudobulk_counts.csv")
        lc, _ = C.cluster_expression(n)
        Z = C.zscore_rows(lc.to_numpy(float))
        pos = {g: i for i, g in enumerate(lc.index)}
        bins = C.expression_bins(lc)
        for tf, (targets, signs) in regs.items():
            if tf not in pos:
                continue
            pairs = [(pos[t], s) for t, s in zip(targets, signs) if t in pos]
            if len(pairs) < MIN_TARGETS:
                continue
            idx = np.array([p for p, _ in pairs])
            sg = np.array([s for _, s in pairs], dtype=float)
            score = (Z[idx] * sg[:, None]).mean(axis=0)
            tf_prof = Z[pos[tf]]
            rho = float(C.spearman_rows(score[None, :], tf_prof)[0])
            draws = C.matched_sets(bins, idx, N_RANDOM, rng)
            W = sp.csr_matrix((np.tile(sg / len(idx), N_RANDOM), (np.repeat(np.arange(N_RANDOM), len(idx)),
                                                                   draws.ravel())), shape=(N_RANDOM, Z.shape[0]))
            null = C.spearman_rows(np.asarray(W @ Z), tf_prof)
            eff, p, mu, sd = C.null_effect(rho, null)
            rows.append({"dataset": ds, "chemistry": chem, "tf": tf, "n_targets": len(idx),
                         "rho_tf_vs_targets": rho, "null_mean": mu, "effect_vs_null_sd": eff, "perm_p": p})
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    tfs = C.tf_list()
    if not tfs or C.collectri().empty:
        out.summary(TITLE, "Skipped: TF annotations missing.",
                    [f"Needs tf/TF_names_v_1.01.txt and tf/collectri.tsv in {C.annotations_dir()}."],
                    ["**Not run**: fetch the annotation files first (running_scripts/fetch_annotations.sh), then rerun "
                     "step 3."], [], [])
        return
    out.used("annotations/tf/TF_names_v_1.01.txt", "annotations/tf/collectri.tsv")
    content, enrich, acts = [], [], []
    for ds in C.DATASETS:
        level = C.gene_level(ds)
        universe = set(level.index[level >= MIN_EXPR_CPM])
        N = len(universe)
        tf_ds = set(C.to_dataset_symbols(sorted(tfs), ds).values()) & universe
        regs = regulons(ds, universe)
        C.log(f"  {ds}: {len(tf_ds)} expressed TFs, {len(regs)} regulons with >= {MIN_TARGETS} expressed targets")
        for name, genes in query_groups(out, ds).items():
            g = set(genes) & universe
            if len(g) < MIN_QUERY:
                continue
            k = len(g & tf_ds)
            content.append({"dataset": ds, "gene_set": name, "n_genes": len(g), "n_tfs": k,
                            "expected": len(g) * len(tf_ds) / N, "p": stats.hypergeom.sf(k - 1, N, len(tf_ds), len(g)),
                            "tfs": "|".join(sorted(g & tf_ds)[:40])})
            ps = []
            for tf, (targets, _) in regs.items():
                t = set(targets)
                kk = len(t & g)
                ps.append({"dataset": ds, "gene_set": name, "tf": tf, "n_targets": len(t), "overlap": kk,
                           "fold": kk / (len(t) * len(g) / N), "p": stats.hypergeom.sf(kk - 1, N, len(t), len(g)),
                           "tf_in_set": tf in g, "targets_in_set": "|".join(sorted(t & g)[:25])})
            ps = pd.DataFrame(ps)
            ps["q"] = C.bh(ps.p)
            enrich.append(ps[(ps.q < 0.25) & (ps.overlap >= 3)])
        acts.append(activity(out, ds, regs, rng))
    content = pd.DataFrame(content)
    if not content.empty:
        content["q"] = C.bh(content.p)
    enrich = pd.concat(enrich, ignore_index=True) if enrich else pd.DataFrame()
    act = pd.concat(acts, ignore_index=True)
    acomb = C.combine_chemistries(act, ["dataset", "tf"], labels=("tracks its targets", "anti-tracks its targets"),
                                  carry=("rho_tf_vs_targets", "n_targets"))
    p03 = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_combined.csv"
    tf_age = pd.DataFrame()
    if p03.exists() and not acomb.empty:
        out.used("results/03_age_trends_within_cell_class/age_trends_combined.csv")
        t03 = pd.read_csv(p03, low_memory=False, dtype={"panels": str})
        tf_age = (t03[(t03.tier.fillna("") != "") & t03.gene.isin(set(acomb.tf))]
                  .assign(trend=lambda d: d.cell_class + " " + d.direction.str.replace(" with age", ""))
                  .groupby(["dataset", "gene"]).trend.apply(lambda s: "; ".join(s)).rename("tf_age_trends").reset_index())
        acomb = acomb.merge(tf_age, left_on=["dataset", "tf"], right_on=["dataset", "gene"], how="left").drop(columns="gene")
    if not enrich.empty and not acomb.empty:
        enrich = enrich.merge(acomb[["dataset", "tf", "tier", "direction"]].rename(
            columns={"tier": "activity_tier", "direction": "activity"}), on=["dataset", "tf"], how="left")
    out.write(content, "tf_content", "Per gene set: TFs among its expressed genes vs the expressed genome")
    out.write(enrich, "regulon_enrichment", "Per gene set x TF: CollecTRI targets in the set (q < 0.25, >= 3); "
              "with the TF's regulon-activity tier")
    out.write(act, "regulon_activity_per_stratum", "Per stratum x TF: Spearman of TF expression with its targets' "
              "signed mean Z across clusters, vs matched random target sets")
    out.write(acomb, "regulon_activity_combined", "Per TF: v2 x v3 combined; tier; the TF's own replicated age trends (03)")

    f = []
    if not content.empty:
        c = content[(content.q < 0.05) & (content.n_tfs > content.expected)]
        f.append("**Gene sets rich in TFs** (TFs / expected; q < 0.05): " + ("; ".join(
            f"{r.dataset} {r.gene_set} {r.n_tfs} / {r.expected:.1f} ({', '.join(r.tfs.split('|')[:6])})"
            for r in c.sort_values("q").head(MAX_LISTED * 2).itertuples()) if len(c) else "none") + ".")
    if not acomb.empty:
        for ds in C.DATASETS:
            a = acomb[(acomb.dataset == ds) & (acomb.tier != "") & (acomb.stouffer_z > 0)]
            f.append(f"**{ds}: regulons whose TF tracks its targets across clusters** ({len(a)} of "
                     f"{int((acomb.dataset == ds).sum())} tested; strongest): " + ", ".join(
                         f"{r.tf} ({r.rho_tf_vs_targets_v2:+.2f}/{r.rho_tf_vs_targets_v3:+.2f})"
                         for r in a.sort_values("stouffer_z", ascending=False).head(12).itertuples()) + ".")
    if not enrich.empty:
        e = enrich[(enrich.q < 0.05) & (enrich.get("activity_tier", "").fillna("") != "")
                   & (enrich.get("activity", "") == "tracks its targets")]
        e = e[e.gene_set.str.startswith(("list:", "seed:", "module:"))]
        parts = []
        for (ds, gs), g in e.groupby(["dataset", "gene_set"], sort=False):
            g = g.sort_values("p")
            parts.append(f"{ds} {gs}: " + ", ".join(f"{r.tf} ({r.overlap}, {r.fold:.1f}x)" for r in g.head(5).itertuples()))
        f.append("**Candidate drivers** (TF targets enriched in the set, q < 0.05, and the TF tracks its targets in "
                 "both donor sets; overlap, fold): " + (" | ".join(parts[:MAX_LISTED * 2]) if parts else "none") + ".")
    out.summary(
        TITLE,
        "Which TFs sit in the modules and NDD lists, whose targets are over-represented in them, and is that TF's "
        "expression tied to its targets' in these data?",
        ["TFs from Lambert et al. 2018; targets from CollecTRI (OmniPath), signed where known; symbols mapped to "
         f"each dataset; universe = expressed genes (mean CPM >= {MIN_EXPR_CPM:g}).",
         f"A: hypergeometric TF content. B: hypergeometric target enrichment for regulons with >= {MIN_TARGETS} "
         "expressed targets, BH within each set.",
         f"C: Spearman of TF expression with its targets' signed mean Z across clusters vs {N_RANDOM} random target "
         "sets matched on mean level x spread; v2 x v3 signed Stouffer, BH, tiered."],
        f,
        ["CollecTRI targets come mostly from other tissues and cell lines; a regulon enriched here is a "
         "hypothesis about regulation in developing brain, not evidence of it.",
         "A TF tracking its targets across clusters can reflect shared cell identity (both high in one cell "
         "type) rather than regulation; within-class tests would be stricter.",
         "TF mRNA is a weak proxy for TF activity."],
        ["Within-class regulon activity (as 07's within_class context); motif or ATAC evidence for direct binding."])


if __name__ == "__main__":
    main()
