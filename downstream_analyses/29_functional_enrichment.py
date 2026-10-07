#!/usr/bin/env python3
"""29 - What do the gene groups found here do? (GO, Reactome, KEGG enrichment)

Question: steps 07-09 found co-expression modules, list sub-modules, genes
that change with age within cell types, and proliferation / phase classes of
genes. Which biological processes, components, functions and pathways do they
carry -- and which do the user lists themselves carry?

Method
  Gene-set libraries from <annotations>/gmt/ (fetch_annotations.sh: GO BP / CC /
  MF 2023, Reactome 2022, KEGG 2021 via Enrichr), term genes mapped to each
  dataset's symbols (exact, case, Ensembl id, the other file, HGNC previous
  symbols). Universe = the dataset's expressed genes (mean CPM >= 1 over cell
  classes and chemistries) -- an enrichment test is only fair against genes the
  data could have found. Terms with 10-500 universe genes. Per query group:
  hypergeometric p for every term, BH within the query; top terms reported.
  Query groups: robust 08 modules; 07 sub-modules; 03 replicated age-trend genes
  per class and direction (>= 15 genes); 09 proliferation classes and S / G2M
  leans; user lists.

Inputs: <annotations>/gmt/*.gmt; results/03, 07, 08, 09 tables; gene lists;
csv_exports/<ds>__<chem>/09_pseudobulk/cell_class__pseudobulk_counts.csv (universe)
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

SLUG = "29_functional_enrichment"
TITLE = "Functional enrichment of modules, sub-modules, age-trend genes and phase classes"
MIN_TERM, MAX_TERM = 10, 500
MIN_QUERY = 10
MIN_EXPR_CPM = 1.0
TOP_TERMS = 5
KEEP_Q = 0.25


def queries(out: C.Output) -> list[tuple[str, str, str, list[str]]]:
    """(dataset, kind, name, genes)."""
    q = []
    p = C.RESULTS / "08_coexpression_modules" / "modules.csv"
    if p.exists():
        out.used("results/08_coexpression_modules/modules.csv")
        m = pd.read_csv(p)
        for r in m[m.robust].itertuples():
            ann = getattr(r, "annotation", "")
            q.append((r.dataset, "module", f"{r.module}" + (f" ({ann})" if isinstance(ann, str) and ann else ""),
                      r.genes.split("|")))
    p = C.RESULTS / "07_gene_list_coherence" / "submodules.csv"
    if p.exists():
        out.used("results/07_gene_list_coherence/submodules.csv")
        s = pd.read_csv(p)
        for r in s[~s.no_split.astype(bool)].itertuples():
            q.append((r.dataset, "sub-module", r.submodule, r.genes.split("|")))
    p = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_combined.csv"
    if p.exists():
        out.used("results/03_age_trends_within_cell_class/age_trends_combined.csv")
        t = pd.read_csv(p, low_memory=False, dtype={"panels": str})
        t = t[t.tier == "replicated"]
        for (ds, cls, d), g in t.groupby(["dataset", "cell_class", "direction"]):
            q.append((ds, "age trend", f"{cls}: {d}", list(g.gene)))
    p = C.RESULTS / "09_cell_cycle_programs" / "gene_phase_map.csv"
    if p.exists():
        out.used("results/09_cell_cycle_programs/gene_phase_map.csv")
        g = pd.read_csv(p)
        for ds, h in g.groupby("dataset"):
            for cls in ("proliferative", "anti-proliferative"):
                q.append((ds, "proliferation", cls, list(h.gene[h.proliferation_class == cls])))
            for lean in ("S", "G2/M"):
                q.append((ds, "proliferation", f"proliferative, {lean}-leaning", list(h.gene[h.phase_lean == lean])))
    if C.gene_lists():
        for ds in C.DATASETS:
            for name, genes in C.mapped_lists(ds).items():
                q.append((ds, "gene list", name, genes))
    return q


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    libs = C.gmt_libraries()
    if not libs:
        out.summary(TITLE, "Skipped: no gene-set libraries.",
                    [f"Needs .gmt files in {C.annotations_dir() / 'gmt'}."],
                    ["**Not run**: fetch the annotation files first (running_scripts/fetch_annotations.sh, on the "
                     "login node), then rerun step 3."], [], [])
        return
    out.used(*[f"annotations/gmt/{k}.gmt" for k in libs])
    qs = queries(out)
    rows = []
    for ds in C.DATASETS:
        level = C.gene_level(ds)
        for chem in C.CHEMISTRIES:
            out.used(f"{C.ns(ds, chem)}/09_pseudobulk/cell_class__pseudobulk_counts.csv")
        universe = list(level.index[level >= MIN_EXPR_CPM])
        upos = {g: i for i, g in enumerate(universe)}
        N = len(universe)
        terms, cols, rws = [], [], []
        for lib, gs in libs.items():
            allg = sorted({g for v in gs.values() for g in v})
            mp = C.to_dataset_symbols(allg, ds)
            for term, genes in gs.items():
                idx = sorted({upos[mp[g]] for g in genes if g in mp and mp[g] in upos})
                if MIN_TERM <= len(idx) <= MAX_TERM:
                    k = len(terms)
                    terms.append((lib, term, len(idx)))
                    rws += [k] * len(idx)
                    cols += idx
        T = sp.csr_matrix((np.ones(len(rws)), (rws, cols)), shape=(len(terms), N))
        tsize = np.array([t[2] for t in terms])
        C.log(f"  {ds}: universe {N:,} genes, {len(terms):,} terms of {MIN_TERM}-{MAX_TERM} genes")
        for d, kind, name, genes in qs:
            if d != ds:
                continue
            qi = sorted({upos[g] for g in genes if g in upos})
            if len(qi) < MIN_QUERY:
                continue
            v = np.zeros(N)
            v[qi] = 1
            k = np.asarray(T @ v).ravel()
            n = len(qi)
            p = stats.hypergeom.sf(k - 1, N, tsize, n)
            qv = C.bh(p)
            exp = tsize * n / N
            keep = (qv < KEEP_Q) & (k >= 3)
            for j in np.nonzero(keep)[0]:
                rows.append({"dataset": ds, "kind": kind, "query": name, "query_genes": n, "library": terms[j][0],
                             "term": terms[j][1], "term_genes": int(tsize[j]), "overlap": int(k[j]),
                             "fold": k[j] / exp[j], "p": p[j], "q": qv[j]})
    res = pd.DataFrame(rows)
    if not res.empty:
        res = res.sort_values(["dataset", "kind", "query", "p"])
    out.write(res, "enrichment", f"Per query group x term: overlap, fold, hypergeometric p, BH q within the query "
              f"(rows with q < {KEEP_Q} and >= 3 genes)")
    top = (res[res.q < 0.05].groupby(["dataset", "kind", "query"], sort=False).head(TOP_TERMS)
           if not res.empty else res)
    out.write(top, "top_terms", f"Top {TOP_TERMS} terms (q < 0.05) per query group")

    f = ["**Libraries**: " + ", ".join(f"{k} ({len(v):,} terms)" for k, v in libs.items())
         + f"; query groups tested: {sum(1 for _ in qs)}."]
    for kind in ("module", "sub-module", "age trend", "proliferation", "gene list"):
        g = top[top.kind == kind] if len(top) else top
        if g.empty:
            continue
        parts = []
        for (ds, qn), h in g.groupby(["dataset", "query"], sort=False):
            parts.append(f"{ds} {qn}: " + "; ".join(
                f"{r.term.split(' (GO:')[0]} ({r.overlap}, {r.fold:.1f}x, q {C.fmt_p(r.q)})" for r in h.head(3).itertuples()))
        f.append(f"**{kind.capitalize()}s -- top terms** (overlap, fold, q): " + " | ".join(parts[:16])
                 + ("" if len(parts) <= 16 else f" | (+{len(parts) - 16} more groups)") + ".")
    out.summary(
        TITLE,
        "Which biological processes, components, functions and pathways do the modules, list sub-modules, "
        "age-trend genes, proliferation classes and user lists carry?",
        ["Libraries: GO BP / CC / MF 2023, Reactome 2022, KEGG 2021 (Enrichr .gmt, fetch_annotations.sh); term "
         "genes mapped to each dataset's symbols.",
         f"Universe: expressed genes (mean CPM >= {MIN_EXPR_CPM:g}); terms of {MIN_TERM}-{MAX_TERM} universe genes; "
         "hypergeometric p, BH within each query group."],
        f,
        ["Terms overlap heavily (GO's hierarchy), so several top terms often describe one signal.",
         "Enrichment ignores expression level within the universe; highly expressed housekeeping processes "
         "(ribosome, translation) can surface for groups of highly expressed genes."],
        ["Collapse redundant GO terms (e.g. by semantic similarity) for a shorter summary."])


if __name__ == "__main__":
    main()
