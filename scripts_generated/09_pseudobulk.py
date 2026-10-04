#!/usr/bin/env python3
"""09 - Pseudobulk expression matrices (the one script that reads X).

This is the bridge from "too big to analyse" to "a CSV you can work with".
It streams the count matrix in row-chunks and, in a *single pass*, fills one
aggregator per grouping, so cell class, cluster, sample and cell-class x age
pseudobulks all come out of the same read of the disk.

Three complementary matrices are written per grouping, because each answers a
different question and the wrong one quietly misleads:
  * pseudobulk_counts  -- summed raw counts; the correct input for a proper
                          negative-binomial DE test (edgeR / DESeq2)
  * mean_lognorm       -- mean log1p(CP10K); comparable across groups whose
                          sequencing depth differs
  * detection_fraction -- share of cells expressing the gene at all, which
                          separates "high in a few cells" from "broadly on"

Genes are capped (default 12k, ranked by total UMIs) to keep CSVs tractable;
every gene in any panel is force-included regardless of rank.

Outputs (csv_exports/<dataset>/09_pseudobulk/)
  <grouping>__pseudobulk_counts.csv     genes x groups, summed raw counts
  <grouping>__mean_lognorm.csv          genes x groups, mean log1p(CP10K)
  <grouping>__detection_fraction.csv    genes x groups, fraction expressing
  <grouping>__cpm.csv                   genes x groups, pseudobulk CPM
  <grouping>__group_summary.csv         cells, depth and genes detected per group
  gene_selection.csv                    which genes were kept and why
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.aggregate import (GroupAggregator, combine_keys, drop_empty_levels,
                           group_codes)
from lib.io_utils import (Manifest, XReader, add_derived_obs_columns,
                          chemistry_mask, exclusion_mask, gene_frame, log,
                          read_obs, read_var,
                          resolve_cluster_columns, resolve_role)
from lib.panels import all_panel_genes

SCRIPT = "09_pseudobulk"
SUBDIR = "09_pseudobulk"


def build_groupings(obs: pd.DataFrame) -> dict[str, pd.Series]:
    """Which partitions of the cells to pseudobulk, resolved per dataset."""
    out: dict[str, pd.Series] = {}
    cls = resolve_role(obs, "cell_class")
    if cls is not None:
        out["cell_class"] = obs[cls].astype("object")
        if "age_pcw" in obs.columns:
            out["cell_class_x_age"] = combine_keys(
                obs.assign(_a=obs["age_pcw"].astype(str)), [cls, "_a"])
        reg = resolve_role(obs, "region")
        if reg is not None and obs[reg].nunique() > 1:
            out["cell_class_x_region"] = combine_keys(obs, [cls, reg])
    for role in ("sample", "donor", "region", "subregion"):
        col = resolve_role(obs, role)
        if col is not None and obs[col].nunique(dropna=True) > 1:
            out[role] = obs[col].astype("object")
    if "age_pcw" in obs.columns:
        out["age"] = obs["age_pcw"].astype(str)
    for col in resolve_cluster_columns(obs):
        out[f"cluster_{col}"] = obs[col].astype("object")
    return out


def select_genes(path, var: pd.DataFrame, genes: pd.DataFrame, top_n: int,
                 chunk_size: int, limit_cells: int | None):
    """Pick which genes to export: top-N by total UMIs, plus all panel genes.

    Uses var['GeneTotalUMIs'] when the file provides it (human_dev), otherwise
    pays for one extra streaming pass to compute per-gene totals (cortex).
    """
    n_genes = len(genes)
    if "GeneTotalUMIs" in var.columns:
        totals = pd.to_numeric(var["GeneTotalUMIs"], errors="coerce").fillna(0).to_numpy()
        log("  gene totals taken from var['GeneTotalUMIs'] (no extra pass needed)")
    else:
        log("  var has no GeneTotalUMIs -- one extra pass to compute gene totals...")
        totals = np.zeros(n_genes, dtype=np.float64)
        with XReader(path) as xr:
            stop_at = min(limit_cells or xr.n_obs, xr.n_obs)
            for start, stop, chunk in xr.iter_chunks(chunk_size):
                if start >= stop_at:
                    break
                totals += np.asarray(chunk.sum(axis=0)).ravel()
        log("  gene totals computed")

    sel = pd.DataFrame({
        "var_index": genes.index.astype(str),
        "symbol": genes["symbol"].astype(str),
        "total_umis": totals,
    })
    sel["rank_by_umis"] = sel["total_umis"].rank(ascending=False, method="first").astype(int)
    panel = all_panel_genes()
    sel["in_panel"] = sel["symbol"].isin(panel)
    sel["selected"] = (sel["rank_by_umis"] <= top_n) | sel["in_panel"]
    # A panel gene with zero counts cannot contribute anything but a zero row.
    sel.loc[sel["total_umis"] <= 0, "selected"] = False
    sel["reason"] = np.where(
        sel["rank_by_umis"] <= top_n,
        np.where(sel["in_panel"], "top_umis+panel", "top_umis"),
        np.where(sel["in_panel"], "panel_rescue", "not_selected"))
    return sel


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    # obs stays FULL length here: X chunks are indexed against every cell, so
    # the chemistry filter is applied to the group codes below (code -1, which
    # GroupAggregator already excludes) rather than by subsetting obs.
    obs_full = add_derived_obs_columns(read_obs(path), key)
    chem_keep = chemistry_mask(obs_full, chem)
    if chem_keep is None:
        log(f"  no chemistry column in this file -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    chem_keep = chem_keep.to_numpy()
    excl_keep, _ = exclusion_mask(obs_full, key, within=chem_keep)
    chem_keep = chem_keep & excl_keep
    if chem_keep.sum() == 0:
        log(f"  no cells with chemistry={chem}; skipping")
        man.flush()
        return
    obs = obs_full
    var = read_var(path)
    genes = gene_frame(var, key)

    with XReader(path) as xr:
        n_obs_total, n_vars = xr.n_obs, xr.n_vars
        encoding = xr.encoding
    log(f"  X: {n_obs_total:,} x {n_vars:,} ({encoding})")
    if n_vars != len(genes):
        raise SystemExit(f"var has {len(genes)} rows but X has {n_vars} columns")
    n_cells = min(args.limit_cells, n_obs_total) if args.limit_cells else n_obs_total
    if args.limit_cells:
        obs = obs.iloc[:n_cells]
        chem_keep = chem_keep[:n_cells]
        log(f"  LIMITED to first {n_cells:,} cells (smoke test)")
    log(f"  chemistry={chem or 'pooled'}: {int(chem_keep.sum()):,} / {n_cells:,} cells")

    sel = select_genes(path, var, genes, args.top_genes, args.chunk_size,
                       args.limit_cells)
    gene_mask = sel["selected"].to_numpy()
    kept_symbols = sel.loc[gene_mask, "symbol"].tolist()
    n_kept = int(gene_mask.sum())
    log(f"  genes selected: {n_kept:,} / {n_vars:,} "
        f"({int(sel['in_panel'].sum())} panel genes, "
        f"{int((sel['reason'] == 'panel_rescue').sum())} rescued below the UMI cut)")
    man.write(sel, "gene_selection",
              "Per gene: total UMIs, rank, panel membership and whether it was exported",
              subdir=SUBDIR)

    groupings = build_groupings(obs)
    log(f"  groupings: {', '.join(groupings)}")

    aggs, codes_map, names_map, skipped = {}, {}, {}, []
    est_bytes = 0
    for name, series in groupings.items():
        codes, levels = group_codes(series)
        # Excluding the other chemistry by code rather than by slicing X keeps
        # the single streaming pass intact; the aggregator skips code -1.
        codes = np.where(chem_keep, codes, -1)
        if (codes >= 0).sum() == 0:
            log(f"  skipping grouping '{name}': no cells left after the chemistry filter")
            continue
        codes, levels = drop_empty_levels(codes, levels)
        if len(levels) > config.MAX_GROUPS_WIDE:
            log(f"  skipping grouping '{name}': {len(levels):,} levels exceeds "
                f"MAX_GROUPS_WIDE={config.MAX_GROUPS_WIDE}")
            skipped.append({"grouping": name, "n_levels": len(levels),
                            "reason": "exceeds MAX_GROUPS_WIDE"})
            continue
        aggs[name] = GroupAggregator(len(levels), n_kept, config.TARGET_SUM)
        codes_map[name], names_map[name] = codes, levels
        est_bytes += len(levels) * n_kept * 8 * 3
    if not aggs:
        log("  no usable groupings")
        man.flush()
        return
    log(f"  accumulator memory ~{est_bytes / 1e9:.2f} GB across {len(aggs)} groupings")

    # ---- the single streaming pass ----------------------------------------
    log(f"  streaming X in chunks of {args.chunk_size:,}...")
    done = 0
    with XReader(path) as xr:
        for start, stop, chunk in xr.iter_chunks(args.chunk_size):
            if start >= n_cells:
                break
            if stop > n_cells:
                chunk = chunk[: n_cells - start]
                stop = n_cells
            totals = np.asarray(chunk.sum(axis=1)).ravel()
            for name, agg in aggs.items():
                agg.update(chunk, codes_map[name][start:stop],
                           totals_chunk=totals, gene_mask=gene_mask)
            done = stop
            log(f"    {done:,} / {n_cells:,} cells ({100 * done / n_cells:.1f}%)")
    log("  streaming complete")

    for name, agg in aggs.items():
        levels = names_map[name]
        man.write(agg.pseudobulk_counts(levels, kept_symbols).reset_index(),
                  f"{name}__pseudobulk_counts",
                  f"Summed raw counts per gene per {name} -- input for edgeR/DESeq2",
                  subdir=SUBDIR)
        man.write(agg.mean_lognorm(levels, kept_symbols).reset_index(),
                  f"{name}__mean_lognorm",
                  f"Mean log1p(CP10K) per gene per {name} -- depth-comparable",
                  subdir=SUBDIR)
        man.write(agg.detection_fraction(levels, kept_symbols).reset_index(),
                  f"{name}__detection_fraction",
                  f"Fraction of cells in each {name} expressing each gene",
                  subdir=SUBDIR)
        man.write(agg.cpm(levels, kept_symbols).reset_index(), f"{name}__cpm",
                  f"Pseudobulk CPM per gene per {name}", subdir=SUBDIR)
        gs = agg.group_summary(levels)
        gs["below_min_cells"] = gs["n_cells"] < config.MIN_CELLS_PER_GROUP
        gs.insert(0, "grouping", name)
        man.write(gs, f"{name}__group_summary",
                  f"Cells, sequencing depth and genes detected per {name}",
                  subdir=SUBDIR)

    if skipped:
        man.write(pd.DataFrame(skipped), "skipped_groupings",
                  "Groupings that exceeded the width cap and were not exported",
                  subdir=SUBDIR)
    man.flush()


def main() -> None:
    parser = cli.build_parser(__doc__)
    parser.add_argument("--top-genes", type=int, default=config.TOP_GENES_PSEUDOBULK,
                        help="max genes to export, ranked by total UMIs "
                             "(panel genes are always added)")
    args = parser.parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
