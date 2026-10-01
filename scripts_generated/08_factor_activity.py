#!/usr/bin/env python3
"""08 - Activity of every latent space per biological group.

Complements script 07: 07 says *what genes* a factor contains, this says
*which cells* use it.  Runs over every multi-dimensional obsm matrix, so it
profiles human_dev's Factors and also cortex's X_scVI / PCA latent spaces.

Two diagnostics matter most:
  * tau specificity across cell classes -- is a factor a cell-identity
    programme or a diffuse one shared by everything?
  * correlation with QC metrics -- a factor tracking total_UMIs or
    fraction_mitochondrial is a technical axis, not biology, and must not be
    interpreted as a module.  Cell-cycle and unspliced-fraction scores are
    also correlated, but flagged separately: proliferation and maturation are
    real biology in a developing brain, not artefacts.  Depth itself differs
    by cell type (neurons carry more UMIs than progenitors), so the technical
    flag is a prompt to check, not a verdict.

Outputs (csv_exports/<dataset>/08_factor_activity/)
  factor_activity_<key>_by_<role>.csv  mean/median activity per group
  factor_specificity_<key>.csv         tau across cell classes + top group
  factor_correlation_<key>.csv         factor x factor correlation
  factor_qc_correlation_<key>.csv      factor x QC correlation + technical flag
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import scipy.sparse as sp

import config
from lib import cli
from lib.io_utils import (Manifest, list_h5ad_keys, load_obs,
                          log, read_elem_at, resolve_cluster_columns,
                          resolve_qc_frame, resolve_role)
from lib.stats_utils import corr_matrix, tau_specificity

SCRIPT = "08_factor_activity"
SUBDIR = "08_factor_activity"
GROUP_ROLES = ["cell_class", "region", "subregion", "chemistry", "donor",
               "sample", "cyclephase"]
MAX_LEVELS = 1000
MAX_DIMS = 200
CORR_SUBSAMPLE = 200_000
TECHNICAL_R = 0.3
# Only these QC roles mark an axis as technical; the rest are biological.
TECHNICAL_ROLES = ["frac_mito", "n_genes", "total_umis", "total_rna",
                   "doublet_score"]
CYCLE_ROLES = ["cell_cycle_score", "cycling_score", "cc_g1", "cc_s", "cc_g2m"]


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    obs, keep = load_obs(path, key, chem, args.limit_cells)
    if obs is None:
        log(f"  no chemistry column in this file -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    if len(obs) == 0:
        log(f"  no cells with chemistry={chem}; skipping")
        man.flush()
        return
    log(f"  {len(obs):,} cells after chemistry filter ({chem or 'pooled'})")
    qc = resolve_qc_frame(obs)
    keys = list_h5ad_keys(path)
    if not keys["obsm"]:
        log("  no obsm matrices present")
        man.flush()
        return

    groupings = {r: resolve_role(obs, r) for r in GROUP_ROLES}
    groupings = {r: c for r, c in groupings.items()
                 if c is not None and obs[c].nunique(dropna=True) <= MAX_LEVELS}
    if "cyclephase_h" in obs.columns:
        groupings["cyclephase"] = "cyclephase_h"
    if "age_pcw" in obs.columns:
        groupings["age"] = "age_pcw"
    for c in resolve_cluster_columns(obs):
        if obs[c].nunique(dropna=True) <= MAX_LEVELS:
            groupings[f"cluster_{c}"] = c
    class_col = resolve_role(obs, "cell_class")
    rng = np.random.default_rng(config.RANDOM_SEED)

    for ekey in keys["obsm"]:
        mat = read_elem_at(path, f"obsm/{ekey}")
        if mat is None:
            continue
        if sp.issparse(mat):
            log(f"  skipping obsm/{ekey}: sparse ({mat.shape}) -- this is a "
                "neighbourhood/indicator matrix (e.g. Milo nhoods), not an embedding")
            continue
        mat = np.asarray(mat)
        if mat.ndim != 2 or mat.shape[1] < 2:
            continue
        if args.limit_cells:
            mat = mat[:args.limit_cells]
        # obsm rows align with the unfiltered cells; apply the same mask.
        if keep is not None and mat.shape[0] == keep.size:
            mat = mat[keep]
        if mat.shape[0] != len(obs):
            log(f"  skipping obsm/{ekey}: {mat.shape[0]} rows != {len(obs)} cells")
            continue
        if mat.shape[1] > MAX_DIMS:
            log(f"  skipping obsm/{ekey}: {mat.shape[1]} dims exceeds cap {MAX_DIMS}")
            continue
        mat = np.nan_to_num(mat.astype(np.float64), nan=0.0)
        dims = [f"{ekey}_{i+1}" for i in range(mat.shape[1])]
        log(f"  {ekey}: {mat.shape[0]:,} cells x {mat.shape[1]} dims")
        F = pd.DataFrame(mat, columns=dims, index=obs.index)

        for role, col in groupings.items():
            lab = obs[col].astype("object")
            g = F.groupby(lab.to_numpy(), observed=True)
            mean = g.mean()
            mean.index.name = "group_level"
            out = mean.reset_index()
            out.insert(1, "n_cells", g.size().reindex(mean.index).to_numpy())
            out.insert(0, "statistic", "mean")
            med = g.median()
            med.index.name = "group_level"
            out_med = med.reset_index()
            out_med.insert(1, "n_cells", g.size().reindex(med.index).to_numpy())
            out_med.insert(0, "statistic", "median")
            man.write(pd.concat([out, out_med], ignore_index=True),
                      f"factor_activity_{ekey}_by_{role}",
                      f"Mean and median {ekey} activity per {col}", subdir=SUBDIR)

        # -- identity programme or diffuse? ---------------------------------
        if class_col is not None:
            by_class = F.groupby(obs[class_col].astype("object").to_numpy(),
                                 observed=True).mean()
            # tau needs non-negative input; shift each factor to its own minimum
            shifted = by_class.to_numpy().T
            shifted = shifted - np.nanmin(shifted, axis=1, keepdims=True)
            spec = pd.DataFrame({
                "factor": dims,
                "tau_across_cell_class": tau_specificity(shifted),
                "top_cell_class": by_class.idxmax().to_numpy(),
                "top_value": by_class.max().to_numpy(),
                "range_across_cell_class": (by_class.max() - by_class.min()).to_numpy(),
            })
            man.write(spec.sort_values("tau_across_cell_class", ascending=False),
                      f"factor_specificity_{ekey}",
                      f"tau specificity of each {ekey} dimension across cell classes "
                      "(1 = restricted to one class, 0 = shared by all)",
                      subdir=SUBDIR)

        # -- correlation structure (subsampled for the big file) -------------
        n = len(F)
        take = (rng.choice(n, CORR_SUBSAMPLE, replace=False) if n > CORR_SUBSAMPLE
                else np.arange(n))
        take.sort()
        cm = corr_matrix(F.iloc[take])
        if not cm.empty:
            cm.index.name = "factor"
            man.write(cm.reset_index(), f"factor_correlation_{ekey}",
                      f"Pearson correlation between {ekey} dimensions across cells "
                      f"(n={len(take):,})", subdir=SUBDIR)

        if not qc.empty:
            joint = pd.concat([F.iloc[take].reset_index(drop=True),
                               qc.iloc[take].reset_index(drop=True)], axis=1)
            full = corr_matrix(joint)
            qc_cols = [c for c in qc.columns if c in full.columns]
            fac_cols = [c for c in dims if c in full.columns]
            if qc_cols and fac_cols:
                block = full.loc[fac_cols, qc_cols].copy()
                block["max_abs_qc_correlation"] = block[qc_cols].abs().max(axis=1)
                tech = [c for c in qc_cols if c in TECHNICAL_ROLES]
                cyc = [c for c in qc_cols if c in CYCLE_ROLES]
                block["max_abs_technical_correlation"] = (
                    block[tech].abs().max(axis=1) if tech else np.nan)
                block["likely_technical"] = block["max_abs_technical_correlation"] > TECHNICAL_R
                block["tracks_cell_cycle"] = (
                    block[cyc].abs().max(axis=1) > TECHNICAL_R if cyc else False)
                block.index.name = "factor"
                man.write(block.reset_index(), f"factor_qc_correlation_{ekey}",
                          f"Correlation of each {ekey} dimension with QC metrics; "
                          f"likely_technical = |r|>{TECHNICAL_R} with depth/mito/"
                          "doublet metrics; tracks_cell_cycle flagged separately "
                          "(biology, not artefact)",
                          subdir=SUBDIR)
                n_tech = int(block["likely_technical"].sum())
                if n_tech:
                    log(f"    {n_tech}/{len(fac_cols)} {ekey} dims flagged as likely technical")
                n_cyc = int(block["tracks_cell_cycle"].sum())
                if n_cyc:
                    log(f"    {n_cyc}/{len(fac_cols)} {ekey} dims track the cell cycle")
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
