#!/usr/bin/env python3
"""13 - Keys for comparing the two datasets.

The cortex file is indexed by gene SYMBOL with unversioned accessions; the
whole-brain file is indexed by VERSIONED Ensembl accession and carries the real
symbols in var['Gene'] (its column literally named `gene_symbol` holds
accessions -- a trap this pipeline corrects for).  Joining them naively on the
var index produces an empty intersection.

This writes the mapping tables that make a cortex-vs-whole-brain comparison
valid: a shared gene list keyed on version-stripped accession, and label
mappings for the categorical columns whose vocabularies differ ('PostM' vs
'Post-M', 5 cell classes vs 12).

Outputs (csv_exports/_cross_dataset/)
  gene_id_map.csv          per dataset: var index, symbol, accession, stripped key
  shared_genes.csv         genes present in both, with both datasets' local ids
  gene_id_overlap_summary.csv  how well each join key performs
  label_overlap_<role>.csv levels of each categorical role per dataset
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import (Manifest, add_derived_obs_columns, gene_frame, log,
                          read_obs, read_var, resolve_role)

SCRIPT = "13_cross_dataset_keys"
CROSS_KEY = "_cross_dataset"
LABEL_ROLES = ["cell_class", "region", "subregion", "chemistry", "cyclephase", "donor"]


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    keys = cli.selected_datasets(args)
    if len(keys) < 2:
        log("cross-dataset comparison needs both datasets; run with --dataset all")
        return

    if getattr(args, "chemistry", "all") != "all":
        log("  (gene-id and label maps are chemistry-independent; running pooled)")
    man = Manifest(CROSS_KEY, SCRIPT)
    frames, labels = {}, {}
    for key in keys:
        path = cli.resolve_h5ad(key)
        log(f"reading {key}...")
        g = gene_frame(read_var(path), key).reset_index()
        frames[key] = g
        obs = add_derived_obs_columns(read_obs(path), key)
        labels[key] = {}
        for role in LABEL_ROLES:
            col = ("cyclephase_h" if role == "cyclephase" and "cyclephase_h" in obs.columns
                   else resolve_role(obs, role))
            if col is not None:
                labels[key][role] = (obs[col].astype(str).value_counts()
                                     .rename_axis("level").reset_index(name="n_cells"))
        log(f"  {key}: {len(g):,} genes, roles found: {', '.join(labels[key])}")

    allg = pd.concat(frames.values(), ignore_index=True)
    man.write(allg, "gene_id_map",
              "Per dataset gene identity: var index, symbol, accession and the "
              "version-stripped accession used as the join key")

    a, b = keys[0], keys[1]
    ga, gb = frames[a], frames[b]
    overlap_rows = []
    for join_key in ("accession_base", "symbol"):
        if join_key not in ga.columns or join_key not in gb.columns:
            continue
        sa = set(ga[join_key].dropna().astype(str)) - {"nan", "<NA>", ""}
        sb = set(gb[join_key].dropna().astype(str)) - {"nan", "<NA>", ""}
        shared = sa & sb
        overlap_rows.append({
            "join_key": join_key,
            f"n_{a}": len(sa), f"n_{b}": len(sb),
            "n_shared": len(shared),
            f"fraction_of_{a}": len(shared) / len(sa) if sa else np.nan,
            f"fraction_of_{b}": len(shared) / len(sb) if sb else np.nan,
        })
        log(f"  join on {join_key}: {len(shared):,} shared genes")
    if overlap_rows:
        man.write(pd.DataFrame(overlap_rows), "gene_id_overlap_summary",
                  "Shared-gene counts under each candidate join key -- "
                  "accession_base should outperform raw symbol matching")

    # Build the actual shared-gene table on the best-performing key.
    best = "accession_base" if any(r["join_key"] == "accession_base" and r["n_shared"] > 0
                                   for r in overlap_rows) else "symbol"
    log(f"  building shared gene table on '{best}'")
    la = ga.dropna(subset=[best]).drop_duplicates(best)
    lb = gb.dropna(subset=[best]).drop_duplicates(best)
    shared = la.merge(lb, on=best, suffixes=(f"_{a}", f"_{b}"))
    cols = [best] + [c for c in shared.columns
                     if c.endswith((f"_{a}", f"_{b}")) and "dataset_" not in c]
    man.write(shared[cols], "shared_genes",
              f"Genes present in both datasets, joined on {best}, with each "
              "dataset's local identifiers for subsetting")

    for role in LABEL_ROLES:
        present = {k: labels[k][role] for k in keys if role in labels[k]}
        if len(present) < 2:
            continue
        merged = None
        for k, df in present.items():
            df = df.rename(columns={"n_cells": f"n_cells_{k}"})
            merged = df if merged is None else merged.merge(df, on="level", how="outer")
        for k in present:
            merged[f"in_{k}"] = merged[f"n_cells_{k}"].notna()
        merged["in_both"] = merged[[f"in_{k}" for k in present]].all(axis=1)
        man.write(merged.sort_values("in_both", ascending=False), f"label_overlap_{role}",
                  f"Levels of '{role}' in each dataset and whether they are shared")
    man.flush()


if __name__ == "__main__":
    main()
