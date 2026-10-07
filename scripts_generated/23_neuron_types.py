#!/usr/bin/env python3
"""23 - Excitatory vs inhibitory identity of every neuron and neuroblast.

The "Neuron" and "Neuroblast" labels pool glutamatergic (excitatory) and
GABAergic / glycinergic (inhibitory) cells, which differ in origin (dorsal vs
ventral progenitors), timing, and in which disease genes they express. This
pass scores every neuron and neuroblast for each identity and exports the
tables step 3 needs to analyse the two separately.

Per cell: summed CP10K of the excitatory genes (SLC17A6, SLC17A7, SLC17A8 --
vesicular glutamate transporters -- and NEUROD2, NEUROD6, TBR1) and of the
inhibitory genes (GAD1, GAD2, SLC32A1, SLC6A5 -- GABA / glycine synthesis and
transport -- and DLX1, DLX2, DLX5, DLX6), each as log1p. Sums rather than
means, so one strongly expressed marker suffices (hindbrain glutamatergic
neurons express SLC17A6 but not the cortical NEUROD6 / TBR1). A cell is
excitatory when its excitatory score exceeds the inhibitory one by MARGIN and
reaches MIN_SCORE, inhibitory the other way round, ambiguous when both reach
MIN_SCORE within MARGIN of each other, and unassigned when neither reaches it
(e.g. monoaminergic or very immature cells).

Streams X once (stage 2). Needs 09's gene selection for the pseudobulks.

Outputs (csv_exports/<dataset>__<chem>/23_neuron_types/)
  neuron_type_marker_check.csv                 mean scores and marker detection per type x class
  neuron_type_by_age.csv / _by_donor.csv       counts per class x type per age / donor
  neuron_type_by_region_x_age.csv              counts per class x type per region x age (if regions)
  neuron_type_x_age__{pseudobulk_counts,group_summary}.csv            'Class:type | age'
  neuron_type_x_region_x_age__{pseudobulk_counts,group_summary}.csv   'Class:type | region | age'
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.aggregate import GroupAggregator
from lib.io_utils import (Manifest, XReader, add_derived_obs_columns, chemistry_mask,
                          exclusion_mask, gene_frame, log, read_obs, read_var, resolve_role)

SCRIPT = "23_neuron_types"
SUBDIR = "23_neuron_types"
CLASSES = ["Neuron", "Neuroblast"]
EXC = ["SLC17A6", "SLC17A7", "SLC17A8", "NEUROD2", "NEUROD6", "TBR1"]
INH = ["GAD1", "GAD2", "SLC32A1", "SLC6A5", "DLX1", "DLX2", "DLX5", "DLX6"]
MIN_SCORE = 1.0     # log1p(summed CP10K): about one marker UMI in a cell of <= 5,000 UMIs
MARGIN = 0.5        # log1p units between the two scores
TYPES = ["excitatory", "inhibitory", "ambiguous", "unassigned"]


def call(e: np.ndarray, i: np.ndarray) -> np.ndarray:
    d = e - i
    return np.where((e >= MIN_SCORE) & (d > MARGIN), "excitatory",
                    np.where((i >= MIN_SCORE) & (d < -MARGIN), "inhibitory",
                             np.where((e >= MIN_SCORE) | (i >= MIN_SCORE), "ambiguous", "unassigned")))


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)
    obs = add_derived_obs_columns(read_obs(path), key)
    keep = chemistry_mask(obs, chem)
    if keep is None:
        log(f"  no chemistry column -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    keep = keep.to_numpy() & exclusion_mask(obs, key, within=keep.to_numpy())[0]
    cls_col = resolve_role(obs, "cell_class")
    if cls_col is None:
        log("  no cell-class column; skipping")
        man.flush()
        return
    with XReader(path) as xr:
        n_obs = xr.n_obs
    n_cells = min(args.limit_cells, n_obs) if args.limit_cells else n_obs
    cls = obs[cls_col].astype(str).to_numpy()
    sel = keep & np.isin(cls, CLASSES)
    sel[n_cells:] = False
    log(f"  {int(sel.sum()):,} neurons / neuroblasts (chemistry={chem or 'pooled'})")
    if sel.sum() == 0:
        man.flush()
        return
    genes = gene_frame(read_var(path), key)
    sym = genes["symbol"].astype(str).to_numpy()
    e_idx = [i for i, s in enumerate(sym) if s in EXC]
    i_idx = [i for i, s in enumerate(sym) if s in INH]
    log(f"  excitatory markers: {', '.join(sym[e_idx])} | inhibitory markers: {', '.join(sym[i_idx])}")
    if not e_idx or not i_idx:
        log("  a marker group is missing; skipping")
        man.flush()
        return

    reg_col = resolve_role(obs, "region")
    sel_path = config.CSV_EXPORTS / ns / "09_pseudobulk" / "gene_selection.csv"
    gene_mask = kept = None
    if sel_path.exists():
        gs = pd.read_csv(sel_path)
        if len(gs) == len(sym):
            gene_mask = gs["selected"].astype(bool).to_numpy()
            kept = gs.loc[gene_mask, "symbol"].astype(str).tolist()
    if gene_mask is None:
        log("  09 gene selection not found -- pseudobulks skipped (run 09 first)")

    rows = {"e": [], "i": [], "idx": [], "ed": [], "id": [], "umi": []}
    with XReader(path) as xr:
        for start, stop, chunk in xr.iter_chunks(args.chunk_size):
            if start >= n_cells:
                break
            stop = min(stop, n_cells)
            chunk = chunk[: stop - start]
            m = sel[start:stop]
            if not m.any():
                continue
            tot = np.asarray(chunk.sum(axis=1)).ravel()[m]
            sub = chunk[m]
            scale = np.divide(config.TARGET_SUM, tot, out=np.zeros_like(tot, dtype=float), where=tot > 0)
            E = np.asarray(sub[:, e_idx].todense()) * scale[:, None]
            I = np.asarray(sub[:, i_idx].todense()) * scale[:, None]
            rows["e"].append(np.log1p(E.sum(axis=1)))
            rows["i"].append(np.log1p(I.sum(axis=1)))
            rows["ed"].append(E > 0)
            rows["id"].append(I > 0)
            rows["umi"].append(tot)
            rows["idx"].append(np.nonzero(m)[0] + start)
    idx = np.concatenate(rows["idx"])
    cells = obs.iloc[idx][[c for c in [cls_col, "age_pcw", resolve_role(obs, "donor"), reg_col] if c]].copy()
    cells.columns = ["cell_class"] + list(cells.columns[1:])
    cells["cell_class"] = cells["cell_class"].astype(str)
    if reg_col:
        cells[reg_col] = cells[reg_col].astype(str)
    cells["exc_score"] = np.concatenate(rows["e"])
    cells["inh_score"] = np.concatenate(rows["i"])
    cells["neuron_type"] = call(cells.exc_score.to_numpy(), cells.inh_score.to_numpy())
    cells["umis"] = np.concatenate(rows["umi"])
    ed, idd = np.concatenate(rows["ed"]), np.concatenate(rows["id"])
    log("  called: " + ", ".join(f"{t} {int((cells.neuron_type == t).sum()):,}" for t in TYPES))

    chk = []
    for (c, t), g in cells.groupby(["cell_class", "neuron_type"]):
        pos = np.nonzero((cells.cell_class == c).to_numpy() & (cells.neuron_type == t).to_numpy())[0]
        r = {"cell_class": c, "neuron_type": t, "n_cells": len(g), "mean_exc_score": g.exc_score.mean(),
             "mean_inh_score": g.inh_score.mean(), "median_umis": g.umis.median()}
        r.update({f"det_{s}": float(ed[pos, k].mean()) for k, s in enumerate(sym[e_idx])})
        r.update({f"det_{s}": float(idd[pos, k].mean()) for k, s in enumerate(sym[i_idx])})
        chk.append(r)
    man.write(pd.DataFrame(chk), "neuron_type_marker_check",
              f"Mean scores and marker detection per class x neuron type (margin {MARGIN}, min {MIN_SCORE})",
              subdir=SUBDIR)

    def counts(by: list[str], name: str, desc: str) -> None:
        if not all(b in cells for b in by):
            return
        ct = cells.groupby(by + ["cell_class", "neuron_type"], observed=True).size().unstack("neuron_type", fill_value=0)
        ct = ct.reindex(columns=TYPES, fill_value=0)
        ct.insert(0, "n_cells", ct.sum(axis=1))
        man.write(ct.reset_index().rename(columns={"age_pcw": "age"}), name, desc, subdir=SUBDIR)

    counts(["age_pcw"], "neuron_type_by_age", "Cells per class x neuron type per age")
    don = resolve_role(obs, "donor")
    if don:
        counts([don], "neuron_type_by_donor", "Cells per class x neuron type per donor")
    if reg_col and cells[reg_col].nunique() > 1:
        counts([reg_col, "age_pcw"], "neuron_type_by_region_x_age", "Cells per class x neuron type per region x age")

    if gene_mask is not None and "age_pcw" in cells:
        lab_age = (cells.cell_class + ":" + cells.neuron_type + " | " + cells.age_pcw.astype(str)).to_numpy()
        groupings = {"neuron_type_x_age": lab_age}
        if reg_col and cells[reg_col].nunique() > 1:
            groupings["neuron_type_x_region_x_age"] = (cells.cell_class + ":" + cells.neuron_type + " | "
                                                       + cells[reg_col].astype(str) + " | "
                                                       + cells.age_pcw.astype(str)).to_numpy()
        aggs, codes, levels = {}, {}, {}
        for name, lab in groupings.items():
            assigned = np.isin(cells.neuron_type.to_numpy(), ["excitatory", "inhibitory"])
            lv = sorted(set(lab[assigned]))
            pos = {v: k for k, v in enumerate(lv)}
            full = np.full(len(obs), -1, dtype=np.int64)
            full[idx[assigned]] = [pos[v] for v in lab[assigned]]
            if len(lv) > config.MAX_GROUPS_WIDE:
                log(f"  skipping {name}: {len(lv)} levels")
                continue
            aggs[name] = GroupAggregator(len(lv), int(gene_mask.sum()), config.TARGET_SUM)
            codes[name], levels[name] = full, lv
        with XReader(path) as xr:
            for start, stop, chunk in xr.iter_chunks(args.chunk_size):
                if start >= n_cells:
                    break
                stop = min(stop, n_cells)
                chunk = chunk[: stop - start]
                if not sel[start:stop].any():
                    continue
                totals = np.asarray(chunk.sum(axis=1)).ravel()
                for name, agg in aggs.items():
                    agg.update(chunk, codes[name][start:stop], totals_chunk=totals, gene_mask=gene_mask)
        for name, agg in aggs.items():
            man.write(agg.pseudobulk_counts(levels[name], kept).reset_index(), f"{name}__pseudobulk_counts",
                      f"Summed raw counts per gene per {name} (excitatory and inhibitory cells only)", subdir=SUBDIR)
            gsum = agg.group_summary(levels[name])
            gsum["below_min_cells"] = gsum["n_cells"] < config.MIN_CELLS_PER_GROUP
            gsum.insert(0, "grouping", name)
            man.write(gsum, f"{name}__group_summary", f"Cells and depth per {name}", subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
