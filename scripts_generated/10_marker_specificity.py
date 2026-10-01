#!/usr/bin/env python3
"""10 - Marker genes and specificity, derived from the pseudobulk CSVs.

Reads only csv_exports/ -- no h5ad.  This is deliberate: once script 09 has run,
every expression question becomes a cheap table operation, which is the whole
point of the CSV layer and the pattern the step-3 analyses should follow.

For each grouping it scores every gene on three axes and combines them:
  * tau specificity across groups (1 = restricted to one group)
  * log2 fold change of the top group against the mean of the rest
  * detection difference -- expressed in what share of cells, in-group vs out

A gene can score high on fold change while being present in a handful of cells;
requiring detection as well is what separates a usable marker from an artefact.
For files lacking varm/Loadings (cortex), the co-expression modules here are the
substitute for script 07's factor modules.

Outputs (csv_exports/<dataset>/10_markers/)
  specificity_<grouping>.csv      per gene: tau, top group, fold change, detection
  top_markers_<grouping>.csv      ranked marker genes per group
  gene_coexpression_<grouping>.csv correlation between top variable genes
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import Manifest, log
from lib.stats_utils import entropy_specificity, tau_specificity

SCRIPT = "10_marker_specificity"
SUBDIR = "10_markers"
TOP_MARKERS_PER_GROUP = 50
COEXPRESSION_GENES = 500
MIN_GROUPS = 2


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    gene_col = df.columns[0]
    df = df.set_index(gene_col)
    df.index.name = "gene"
    return df.apply(pd.to_numeric, errors="coerce")


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    man = Manifest(ns, SCRIPT)
    pb_dir = config.CSV_EXPORTS / ns / "09_pseudobulk"
    if not pb_dir.exists():
        log(f"  {pb_dir} not found -- run 09_pseudobulk.py first")
        man.flush()
        return

    matrices = sorted(pb_dir.glob("*__mean_lognorm.csv"))
    if not matrices:
        log("  no *__mean_lognorm.csv found -- run 09_pseudobulk.py first")
        man.flush()
        return
    log(f"  found {len(matrices)} pseudobulk matrices")

    for mpath in matrices:
        grouping = mpath.name.replace("__mean_lognorm.csv", "")
        expr = load_matrix(mpath)
        if expr.shape[1] < MIN_GROUPS:
            log(f"  skipping {grouping}: only {expr.shape[1]} group(s)")
            continue

        dpath = pb_dir / f"{grouping}__detection_fraction.csv"
        detect = load_matrix(dpath) if dpath.exists() else None
        # Small groups make per-group means unstable; drop them from marker
        # calling rather than letting a 3-cell group define a marker.
        spath = pb_dir / f"{grouping}__group_summary.csv"
        if spath.exists():
            gs = pd.read_csv(spath)
            good = gs.loc[gs["n_cells"] >= config.MIN_CELLS_PER_GROUP, "group"].astype(str)
            keep = [c for c in expr.columns if c in set(good)]
            if len(keep) >= MIN_GROUPS and len(keep) < expr.shape[1]:
                log(f"  {grouping}: dropping {expr.shape[1] - len(keep)} group(s) "
                    f"below {config.MIN_CELLS_PER_GROUP} cells")
                expr = expr[keep]
                if detect is not None:
                    detect = detect[[c for c in keep if c in detect.columns]]

        values = expr.to_numpy(dtype=float)
        n_groups = values.shape[1]
        top_idx = np.nanargmax(np.nan_to_num(values, nan=-np.inf), axis=1)
        top_val = values[np.arange(len(values)), top_idx]
        row_sum = np.nansum(values, axis=1)
        rest_mean = (row_sum - np.nan_to_num(top_val)) / max(n_groups - 1, 1)

        spec = pd.DataFrame({
            "gene": expr.index.astype(str),
            "top_group": expr.columns[top_idx],
            "top_mean_lognorm": top_val,
            "rest_mean_lognorm": rest_mean,
            # +eps keeps a silent-elsewhere gene finite instead of inf
            "log2fc_vs_rest": np.log2((top_val + 1e-9) / (rest_mean + 1e-9)),
            "tau_specificity": tau_specificity(values),
            "entropy_specificity": entropy_specificity(values),
            "mean_across_groups": np.nanmean(values, axis=1),
            "sd_across_groups": np.nanstd(values, axis=1),
            "n_groups_expressed": (values > 0).sum(axis=1),
        })
        if detect is not None:
            dv = detect.reindex(index=expr.index, columns=expr.columns).to_numpy(dtype=float)
            d_top = dv[np.arange(len(dv)), top_idx]
            d_rest = (np.nansum(dv, axis=1) - np.nan_to_num(d_top)) / max(n_groups - 1, 1)
            spec["detection_in_top_group"] = d_top
            spec["detection_in_rest"] = d_rest
            spec["detection_difference"] = d_top - d_rest
        spec.insert(0, "grouping", grouping)
        man.write(spec.sort_values("tau_specificity", ascending=False),
                  f"specificity_{grouping}",
                  f"Per-gene specificity across {grouping}: tau, fold change, detection",
                  subdir=SUBDIR)

        # -- ranked markers per group ---------------------------------------
        score = spec["log2fc_vs_rest"].fillna(0) * spec["tau_specificity"].fillna(0)
        if "detection_difference" in spec.columns:
            score = score * spec["detection_difference"].clip(lower=0).fillna(0)
        ranked = spec.assign(marker_score=score)
        tops = (ranked.sort_values("marker_score", ascending=False)
                .groupby("top_group", group_keys=False)
                .head(TOP_MARKERS_PER_GROUP)
                .sort_values(["top_group", "marker_score"], ascending=[True, False]))
        tops["rank_in_group"] = tops.groupby("top_group").cumcount() + 1
        man.write(tops, f"top_markers_{grouping}",
                  f"Top {TOP_MARKERS_PER_GROUP} markers per {grouping} group, ranked by "
                  "specificity x fold change x detection gain", subdir=SUBDIR)

        # -- co-expression modules (the stand-in for varm/Loadings) ---------
        if n_groups >= 4:
            var_rank = expr.var(axis=1).sort_values(ascending=False)
            sub = expr.loc[var_rank.head(min(COEXPRESSION_GENES, len(var_rank))).index]
            sub = sub[sub.var(axis=1) > 0]
            if len(sub) >= 10:
                cm = np.corrcoef(sub.to_numpy())
                cdf = pd.DataFrame(cm, index=sub.index, columns=sub.index)
                cdf.index.name = "gene"
                man.write(cdf.reset_index(), f"gene_coexpression_{grouping}",
                          f"Correlation across {grouping} groups between the "
                          f"{len(sub)} most variable genes -- co-expression modules",
                          subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
