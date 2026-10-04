#!/usr/bin/env python3
"""12 - Spliced / unspliced structure per gene per group (cortex only).

The cortex file carries velocyto layers (spliced, unspliced, ambiguous).  Full
RNA velocity needs the whole matrix, but the quantity velocity is built on --
the unspliced fraction, a proxy for how actively a gene is being transcribed
right now -- aggregates perfectly well into a pseudobulk table.

A gene whose unspliced fraction is high in progenitors and low in neurons is
being switched off; the reverse is being switched on.  That is a directional
signal the plain expression tables cannot give, and it costs one more pass.

Skips cleanly on files without these layers (human_dev has none).

Outputs (csv_exports/<dataset>/12_splicing/)
  <grouping>__spliced_counts.csv     genes x groups, summed spliced counts
  <grouping>__unspliced_counts.csv   genes x groups, summed unspliced counts
  <grouping>__unspliced_ratio.csv    unspliced / (spliced + unspliced)
  unspliced_ratio_summary.csv        per gene: spread of unspliced ratio
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.aggregate import GroupAggregator, drop_empty_levels, group_codes
from lib.io_utils import (Manifest, XReader, add_derived_obs_columns,
                          chemistry_mask, exclusion_mask, gene_frame,
                          list_h5ad_keys, log,
                          read_obs, read_var,
                          resolve_cluster_columns, resolve_role)

SCRIPT = "12_splicing_layers"
SUBDIR = "12_splicing"
SPLICED_KEYS = ["spliced", "Spliced"]
UNSPLICED_KEYS = ["unspliced", "Unspliced"]
MIN_COUNTS_FOR_RATIO = 50


def stream_layer(path, layer: str, codes: np.ndarray, n_groups: int,
                 gene_mask: np.ndarray, n_kept: int, chunk_size: int,
                 n_cells: int) -> GroupAggregator:
    agg = GroupAggregator(n_groups, n_kept, config.TARGET_SUM)
    with XReader(path, key=f"layers/{layer}") as xr:
        for start, stop, chunk in xr.iter_chunks(chunk_size):
            if start >= n_cells:
                break
            if stop > n_cells:
                chunk = chunk[: n_cells - start]
                stop = n_cells
            totals = np.asarray(chunk.sum(axis=1)).ravel()
            agg.update(chunk, codes[start:stop], totals_chunk=totals,
                       gene_mask=gene_mask)
    return agg


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    layers = list_h5ad_keys(path)["layers"]
    spliced = next((k for k in SPLICED_KEYS if k in layers), None)
    unspliced = next((k for k in UNSPLICED_KEYS if k in layers), None)
    if spliced is None or unspliced is None:
        log(f"  no spliced/unspliced layers (layers present: {layers or 'none'}); skipping")
        man.flush()
        return
    log(f"  using layers '{spliced}' and '{unspliced}'")

    # Same reasoning as script 09: keep obs full length and exclude the other
    # chemistry through the group codes, so the layer pass stays one sweep.
    obs = add_derived_obs_columns(read_obs(path), key)
    chem_keep = chemistry_mask(obs, chem)
    if chem_keep is None:
        log(f"  no chemistry column in this file -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    chem_keep = chem_keep.to_numpy()
    excl_keep, _ = exclusion_mask(obs, key, within=chem_keep)
    chem_keep = chem_keep & excl_keep
    if chem_keep.sum() == 0:
        log(f"  no cells with chemistry={chem}; skipping")
        man.flush()
        return
    var = read_var(path)
    genes = gene_frame(var, key)
    with XReader(path) as xr:
        n_obs_total = xr.n_obs
    n_cells = min(args.limit_cells, n_obs_total) if args.limit_cells else n_obs_total
    if args.limit_cells:
        obs = obs.iloc[:n_cells]
        chem_keep = chem_keep[:n_cells]
    log(f"  chemistry={chem or 'pooled'}: {int(chem_keep.sum()):,} / {n_cells:,} cells")

    # Reuse script 09's gene selection so the tables line up gene-for-gene.
    sel_path = config.CSV_EXPORTS / ns / "09_pseudobulk" / "gene_selection.csv"
    if sel_path.exists():
        sel = pd.read_csv(sel_path)
        gene_mask = sel["selected"].to_numpy(dtype=bool)
        log(f"  reusing gene selection from 09_pseudobulk ({gene_mask.sum():,} genes)")
    else:
        log("  09_pseudobulk gene_selection.csv not found; using all genes")
        gene_mask = np.ones(len(genes), dtype=bool)
    kept = genes.loc[gene_mask, "symbol"].astype(str).tolist()
    n_kept = int(gene_mask.sum())

    groupings: dict[str, pd.Series] = {}
    cls = resolve_role(obs, "cell_class")
    if cls is not None:
        groupings["cell_class"] = obs[cls].astype("object")
        if "age_pcw" in obs.columns:
            groupings["cell_class_x_age"] = (obs[cls].astype(str) + " | "
                                             + obs["age_pcw"].astype(str))
    for c in resolve_cluster_columns(obs)[:2]:
        groupings[f"cluster_{c}"] = obs[c].astype("object")
    if not groupings:
        log("  no usable groupings")
        man.flush()
        return

    ratio_rows = []
    for name, series in groupings.items():
        codes, levels = group_codes(series)
        codes = np.where(chem_keep, codes, -1)
        if (codes >= 0).sum() == 0:
            continue
        codes, levels = drop_empty_levels(codes, levels)
        if len(levels) > config.MAX_GROUPS_WIDE:
            continue
        log(f"  {name}: {len(levels)} groups -- streaming both layers...")
        s_agg = stream_layer(path, spliced, codes, len(levels), gene_mask, n_kept,
                             args.chunk_size, n_cells)
        u_agg = stream_layer(path, unspliced, codes, len(levels), gene_mask, n_kept,
                             args.chunk_size, n_cells)

        s_df = s_agg.pseudobulk_counts(levels, kept)
        u_df = u_agg.pseudobulk_counts(levels, kept)
        total = s_df + u_df
        # A ratio from a handful of reads is noise, so mask thin cells outright.
        ratio = (u_df / total.where(total >= MIN_COUNTS_FOR_RATIO)).astype(float)

        man.write(s_df.reset_index(), f"{name}__spliced_counts",
                  f"Summed spliced counts per gene per {name}", subdir=SUBDIR)
        man.write(u_df.reset_index(), f"{name}__unspliced_counts",
                  f"Summed unspliced counts per gene per {name}", subdir=SUBDIR)
        man.write(ratio.reset_index(), f"{name}__unspliced_ratio",
                  f"Unspliced fraction per gene per {name} (NaN below "
                  f"{MIN_COUNTS_FOR_RATIO} counts) -- proxy for active transcription",
                  subdir=SUBDIR)

        r = ratio.to_numpy(dtype=float)
        with np.errstate(invalid="ignore"):
            ratio_rows.append(pd.DataFrame({
                "grouping": name,
                "gene": ratio.index.astype(str),
                "n_groups_measured": np.isfinite(r).sum(axis=1),
                "unspliced_ratio_mean": np.nanmean(r, axis=1),
                "unspliced_ratio_min": np.nanmin(r, axis=1),
                "unspliced_ratio_max": np.nanmax(r, axis=1),
                "unspliced_ratio_range": np.nanmax(r, axis=1) - np.nanmin(r, axis=1),
                "total_counts": total.to_numpy(dtype=float).sum(axis=1),
            }))

    if ratio_rows:
        summary = pd.concat(ratio_rows, ignore_index=True)
        man.write(summary.sort_values("unspliced_ratio_range", ascending=False),
                  "unspliced_ratio_summary",
                  "Per gene per grouping: spread of unspliced fraction across groups -- "
                  "a large range means the gene is being switched on or off",
                  subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
