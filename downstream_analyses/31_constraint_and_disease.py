#!/usr/bin/env python3
"""31 - Mutation intolerance (gnomAD) and Mendelian disease genes in modules and lists.

Question: are the NDD lists, list sub-modules and co-expression modules made
of genes intolerant of loss-of-function mutations (low LOEUF) and of known
Mendelian disease genes -- more than expression-matched genes? And do
constrained genes behave differently in these data: tighter between donors
(15), more connected within their list (07), more often changing with age (03)?

Method
  gnomAD v4.1 constraint (LOEUF = upper bound of the LoF observed/expected
  ratio; MANE / canonical transcript); HPO gene-to-disease (OMIM / Orphanet
  Mendelian genes). Symbols mapped to each dataset's.
  A. Per gene set: mean LOEUF, share in the most constrained decile (LOEUF of
     the genome's bottom 10%) and share of Mendelian disease genes, each against
     2,000 random sets matched on expression decile (constraint and expression
     are correlated).
  B. Per class and stratum: Spearman of LOEUF with 15's variability percentile
     (positive = constrained genes vary less between donors); with 07's
     within-list connectivity; and with |combined Z| of 03's age trends.

Inputs: <annotations>/{constraint/gnomad.v4.1.constraint_metrics.tsv, disease/genes_to_disease.txt};
results/03, 07, 08, 15 tables; gene lists; cell_class pseudobulks (expression level)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "31_constraint_and_disease"
TITLE = "Mutation intolerance (gnomAD LOEUF) and Mendelian disease genes in lists and modules"
N_RANDOM = 2000
MIN_GENES = 5
MAX_LISTED = 10


def per_dataset_gene_table(ds: str, cons: pd.DataFrame, dis: pd.DataFrame) -> pd.DataFrame:
    level = C.gene_level(ds)
    mp = C.to_dataset_symbols(sorted(set(cons.index)), ds)
    c = cons.rename(index=mp)
    c = c[c.index.isin(level.index) & (c.index != "")]
    c = c[~c.index.duplicated()]
    dmp = C.to_dataset_symbols(sorted(set(dis.gene.astype(str))), ds)
    mend = {dmp[g] for g in dis.loc[dis.association_type == "MENDELIAN", "gene"].astype(str) if g in dmp}
    t = pd.DataFrame({"level": np.log2(level + 1)})
    t = t.join(c, how="left")
    t["mendelian"] = t.index.isin(mend).astype(float)
    return t


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    cons, dis = C.gnomad_constraint(), C.disease_genes()
    if cons.empty:
        out.summary(TITLE, "Skipped: gnomAD constraint missing.",
                    [f"Needs constraint/gnomad.v4.1.constraint_metrics.tsv in {C.annotations_dir()}."],
                    ["**Not run**: fetch the annotation files first (running_scripts/fetch_annotations.sh), then rerun "
                     "step 3."], [], [])
        return
    out.used("annotations/constraint/gnomad.v4.1.constraint_metrics.tsv", "annotations/disease/genes_to_disease.txt")
    sets_rows, rel_rows, tables = [], [], {}
    for ds in C.DATASETS:
        for chem in C.CHEMISTRIES:
            out.used(f"{C.ns(ds, chem)}/09_pseudobulk/cell_class__pseudobulk_counts.csv")
        t = per_dataset_gene_table(ds, cons, dis)
        tables[ds] = t
        thr = float(cons.loeuf.quantile(0.10))
        t["top_decile"] = (t.loeuf <= thr).astype(float).where(t.loeuf.notna())
        sets = C.analysis_gene_sets(ds, groups=("ndd",), modules=True, out=out)
        p7 = C.RESULTS / "07_gene_list_coherence" / "submodules.csv"
        if p7.exists():
            sm = pd.read_csv(p7)
            for r in sm[(sm.dataset == ds) & ~sm.no_split.astype(bool)].itertuples():
                sets[f"submodule:{r.submodule}"] = r.genes.split("|")
        has = t.dropna(subset=["loeuf"])
        for stat, label in (("loeuf", "mean LOEUF"), ("top_decile", "share in most constrained decile"),
                            ("mendelian", "share Mendelian disease genes")):
            src = has if stat != "mendelian" else t
            st = C.set_shift_test(src[stat], src.level, sets, rng, N_RANDOM, min_genes=MIN_GENES)
            if not st.empty:
                sets_rows.append(st.assign(dataset=ds, statistic=label))
        # B. relations with step-3 gene properties
        p15 = C.RESULTS / "15_donor_variability" / "gene_variability.csv"
        if p15.exists():
            out.used("results/15_donor_variability/gene_variability.csv")
            v = pd.read_csv(p15)
            for (chem, cls), g in v[v.dataset == ds].groupby(["chemistry", "cell_class"]):
                g = g.set_index("gene").join(has.loeuf, how="inner")
                if len(g) > 100:
                    rel_rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls,
                                     "relation": "LOEUF vs between-donor variability percentile (15)",
                                     "n_genes": len(g), "spearman": float(g.loeuf.corr(g.variability_percentile,
                                                                                        method="spearman"))})
        p07 = C.RESULTS / "07_gene_list_coherence" / "gene_connectivity.csv"
        if p07.exists():
            out.used("results/07_gene_list_coherence/gene_connectivity.csv")
            k = pd.read_csv(p07)
            kc = [c for c in k.columns if c.startswith("connectivity")]
            if kc and "gene" in k:
                k = k[k.dataset == ds] if "dataset" in k else k
                k = k.assign(conn=k[kc].mean(axis=1))          # mean over chemistries
                g = k.groupby("gene").conn.mean().to_frame().join(has.loeuf, how="inner")
                if len(g) > 30:
                    rel_rows.append({"dataset": ds, "chemistry": "both", "cell_class": "all",
                                     "relation": "LOEUF vs within-list connectivity (07)", "n_genes": len(g),
                                     "spearman": float(g.loeuf.corr(g.conn, method="spearman"))})
        p03 = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_combined.csv"
        if p03.exists():
            out.used("results/03_age_trends_within_cell_class/age_trends_combined.csv")
            a = pd.read_csv(p03, low_memory=False, dtype={"panels": str})
            for cls, g in a[a.dataset == ds].groupby("cell_class"):
                g = g.set_index("gene").join(has.loeuf, how="inner")
                if len(g) > 100:
                    rel_rows.append({"dataset": ds, "chemistry": "both", "cell_class": cls,
                                     "relation": "LOEUF vs |age-trend Z| (03)", "n_genes": len(g),
                                     "spearman": float(g.loeuf.corr(g.stouffer_z.abs(), method="spearman"))})
    res = pd.concat(sets_rows, ignore_index=True) if sets_rows else pd.DataFrame()
    if not res.empty:
        res["q"] = np.nan
        for _, ix in res.groupby(["dataset", "statistic"]).groups.items():
            res.loc[ix, "q"] = C.bh(res.loc[ix, "perm_p"])
    rel = pd.DataFrame(rel_rows)
    out.write(res, "set_constraint", "Per gene set x statistic: mean LOEUF / constrained-decile share / Mendelian share "
              "vs expression-matched random sets; BH per dataset x statistic")
    out.write(rel, "constraint_relations", "Spearman of LOEUF with step-3 gene properties (variability, connectivity, "
              "age-trend strength)")
    out.write(pd.concat([t.assign(dataset=ds) for ds, t in tables.items()]).reset_index(names="gene"),
              "gene_constraint", "Per gene: expression level (log2 CPM + 1), LOEUF, pLI, missense Z, Mendelian flag")

    f = []
    if not res.empty:
        for stat in ("mean LOEUF", "share Mendelian disease genes"):
            r = res[(res.statistic == stat) & (res.q < 0.05)]
            low = r[r.effect_vs_null_sd < 0] if stat == "mean LOEUF" else r[r.effect_vs_null_sd > 0]
            label = "more constrained than matched genes (mean LOEUF; null)" if stat == "mean LOEUF" else \
                "richer in Mendelian disease genes (share; null)"
            f.append(f"**Gene sets {label}**: " + ("; ".join(
                f"{x.dataset} {x.gene_set} {x.mean:.2f} vs {x.null_mean:.2f}"
                for x in low.sort_values("effect_vs_null_sd", ascending=stat != "mean LOEUF").head(MAX_LISTED * 2).itertuples())
                if len(low) else "none at q < 0.05") + ".")
    if not rel.empty:
        for relation, g in rel.groupby("relation", sort=False):
            f.append(f"**{relation}** (Spearman; median over classes / strata, range): median {g.spearman.median():+.2f} "
                     f"({g.spearman.min():+.2f} to {g.spearman.max():+.2f}; n = {len(g)}).")
    out.summary(
        TITLE,
        "Are NDD lists, sub-modules and modules enriched for loss-of-function-intolerant and Mendelian disease "
        "genes beyond expression-matched genes, and do constrained genes vary less between donors, connect "
        "more within lists, or change more with age?",
        ["gnomAD v4.1 constraint (LOEUF; MANE / canonical transcript); HPO gene-to-disease (Mendelian = OMIM / "
         "Orphanet genes annotated MENDELIAN); symbols mapped to each dataset.",
         f"A: per set, mean LOEUF, share in the genome's most constrained LOEUF decile, Mendelian share, each vs "
         f"{N_RANDOM:,} random sets matched on expression decile; BH per dataset x statistic.",
         "B: Spearman of LOEUF with 15's variability percentile, 07's connectivity and 03's |combined Z|."],
        f,
        ["Seed panels and many user lists were built from disease genes, so their Mendelian share is circular; "
         "the informative cases are modules and sub-modules found without lists.",
         "LOEUF is unreliable for short genes (few expected LoF variants); its upper bound is conservative there."],
        ["Weight list genes by constraint in 06/07, or test constrained and unconstrained halves separately."])


if __name__ == "__main__":
    main()
