#!/usr/bin/env python3
"""33 - Gene-list entries recovered through HGNC, and what is still missing.

Question: gene lists name genes in whatever symbol version their source used;
the two files use two other versions. How many list entries are matched only
through HGNC previous symbols or aliases (the last route in _common.map_genes,
active once the HGNC table is fetched), and why are the remaining entries
missing -- unknown to HGNC, or known but absent from the file?

Method
  Per list and dataset: entries by match route (exact, case, Ensembl id, the
  other file's annotation, HGNC previous / alias); for unmatched entries,
  whether HGNC knows the symbol (current, previous or alias). Every list
  analysis (06-33) already uses the HGNC route; this report shows its effect.

Inputs: <annotations>/hgnc/hgnc_complete_set.txt; gene lists; csv_exports/_cross_dataset/gene_id_map.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

import _common as C

SLUG = "33_symbol_rescue"
TITLE = "Gene-list symbols recovered through HGNC, and entries still missing"
MAX_EXAMPLES = 12


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    hg = C.hgnc_to_ensembl()
    lists = C.gene_lists()
    if not hg or not lists:
        why = "no HGNC table" if not hg else "no gene lists"
        out.summary(TITLE, f"Skipped: {why}.",
                    [f"Needs hgnc/hgnc_complete_set.txt in {C.annotations_dir()} and gene lists."],
                    [f"**Not run**: {why}. Fetch annotations with running_scripts/fetch_annotations.sh."], [], [])
        return
    out.used("annotations/hgnc/hgnc_complete_set.txt", "_cross_dataset/gene_id_map.csv")
    rows, missing = [], []
    for ds in C.DATASETS:
        for name, genes in lists.items():
            m = C.map_genes(genes, ds)
            counts = m.match.value_counts()
            miss = m[m.match == "missing"]
            known = miss.input.map(lambda x: str(x).upper() in hg or str(x).split(".")[0].startswith("ENSG")).astype(bool)
            rows.append({"dataset": ds, "gene_list": name, "n_entries": len(m),
                         **{f"matched_{k}": int(counts.get(k, 0)) for k in
                            ("exact", "case", "ensembl", "other_file_symbol", "hgnc_previous_or_alias")},
                         "missing_known_to_hgnc": int(known.sum()), "missing_unknown": int((~known).sum()),
                         "rescued_examples": "|".join(f"{a}->{b}" for a, b in
                                                      m.loc[m.match == "hgnc_previous_or_alias", ["input", "symbol"]]
                                                      .head(MAX_EXAMPLES).itertuples(index=False))})
            missing.append(miss.assign(dataset=ds, gene_list=name,
                                       reason=["known to HGNC, not in this file" if k else "not an HGNC symbol"
                                               for k in known]))
    res = pd.DataFrame(rows)
    miss = pd.concat(missing, ignore_index=True) if missing else pd.DataFrame()
    out.write(res, "list_matching_routes", "Per list x dataset: entries matched by each route; unmatched entries "
              "known / unknown to HGNC; examples rescued through HGNC")
    out.write(miss, "unmatched_entries", "Every unmatched list entry and why")
    f = []
    for ds, g in res.groupby("dataset"):
        tot = g.n_entries.sum()
        f.append(f"**{ds}**: {int(g.matched_hgnc_previous_or_alias.sum())} of {tot} list entries matched only through "
                 f"HGNC previous symbols / aliases; still missing {int(g.missing_known_to_hgnc.sum())} known to HGNC "
                 f"but absent from the file, {int(g.missing_unknown.sum())} not HGNC symbols. Per list (rescued / "
                 "missing): " + ", ".join(f"{r.gene_list} {r.matched_hgnc_previous_or_alias} / "
                                          f"{r.missing_known_to_hgnc + r.missing_unknown}" for r in g.itertuples()) + ".")
    ex = res[res.rescued_examples != ""].drop_duplicates("gene_list")
    if len(ex):
        f.append("**Examples rescued**: " + "; ".join(f"{r.gene_list}: {r.rescued_examples.replace('|', ', ')}"
                                                     for r in ex.head(6).itertuples()) + ".")
    out.summary(
        TITLE,
        "How many gene-list entries are matched only through HGNC previous symbols or aliases, and why are the "
        "rest missing?",
        ["_common.map_genes routes: exact symbol, case-insensitive, Ensembl id, the other file's annotation, then "
         "HGNC previous symbols and unambiguous aliases via their Ensembl id (fetch_annotations.sh).",
         "Unmatched entries checked against HGNC current, previous and alias symbols."],
        f,
        ["An alias shared by two genes is never used (ambiguous); such entries stay missing.",
         "'Known to HGNC, not in this file' covers genes the file's annotation omits (often non-coding)."],
        ["Fix symbols at the source lists to the HGNC version the analyses report."])


if __name__ == "__main__":
    main()
