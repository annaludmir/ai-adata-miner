#!/usr/bin/env python3
"""21 - Milo neighbourhoods: cells per neighbourhood x donor, and what each holds.

cortex stores a Milo neighbourhood matrix (obsm['nhoods']: cells x
neighbourhoods, 1 where a cell belongs to a neighbourhood). Neighbourhoods are
small, overlapping groups of transcriptionally similar cells, so abundance can
be tested at finer resolution than clusters and without fixed boundaries.
This pass counts each neighbourhood's cells per donor -- the replicate unit --
and annotates it; step 3 tests abundance against age.

Reads obs and obsm only. Files without an obsm 'nhoods' matrix are skipped.

Outputs (csv_exports/<dataset>__<chem>/21_milo/)
  nhood_counts_by_donor.csv   neighbourhood x donor cell counts
  nhood_annotation.csv        size, dominant class (+ share), dominant cluster of
                              every clustering, mean age, cycling share
  donor_totals.csv            cells per donor x class (denominators)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import scipy.sparse as sp

from lib import cli
from lib.aggregate import group_codes, onehot_csr
from lib.io_utils import (Manifest, list_h5ad_keys, load_obs, log, read_elem_at,
                          resolve_cluster_columns, resolve_role)

SCRIPT = "21_milo_nhoods"
SUBDIR = "21_milo"
NHOOD_KEY = "nhoods"


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)
    if NHOOD_KEY not in list_h5ad_keys(path)["obsm"]:
        log(f"  no obsm/{NHOOD_KEY} -- skipping")
        man.flush()
        return
    obs, keep = load_obs(path, key, chem, args.limit_cells)
    if obs is None or len(obs) == 0:
        log("  no chemistry column or no cells -- skipping")
        man.flush()
        return
    donor_col, cls_col = resolve_role(obs, "donor"), resolve_role(obs, "cell_class")
    if donor_col is None or cls_col is None:
        log("  needs donor and cell-class columns -- skipping")
        man.flush()
        return
    N = read_elem_at(path, f"obsm/{NHOOD_KEY}")
    N = sp.csr_matrix(N) if sp.issparse(N) else sp.csr_matrix(np.asarray(N))
    if args.limit_cells:
        N = N[: args.limit_cells]
    N = N[np.nonzero(keep)[0]]
    N.data = (N.data > 0).astype(np.float64)
    nh = N.shape[1]
    log(f"  {nh:,} neighbourhoods, {len(obs):,} cells")

    def counts_by(col: str) -> pd.DataFrame:
        codes, levels = group_codes(obs[col])
        return pd.DataFrame(np.asarray((N.T @ onehot_csr(codes, len(levels))).todense()),
                            columns=levels)

    by_donor = counts_by(donor_col)
    by_donor.insert(0, "nhood", np.arange(nh))
    size = by_donor.drop(columns="nhood").sum(axis=1)
    keep_nh = size > 0
    man.write(by_donor[keep_nh], "nhood_counts_by_donor", "Cells per neighbourhood x donor", subdir=SUBDIR)

    by_cls = counts_by(cls_col)
    ann = pd.DataFrame({"nhood": np.arange(nh), "n_cells": size.to_numpy()})
    ann["dominant_cell_class"] = by_cls.idxmax(axis=1).to_numpy()
    ann["dominant_class_share"] = (by_cls.max(axis=1) / size.replace(0, np.nan)).to_numpy()
    for col in resolve_cluster_columns(obs):      # every clustering: some are per-chemistry, step 3 picks
        by_cl = counts_by(col)
        ann[f"dominant_{col}"] = by_cl.idxmax(axis=1).to_numpy()
    if "age_pcw" in obs:
        a = obs["age_pcw"].fillna(obs["age_pcw"].median()).to_numpy(float)
        ann["mean_age_pcw"] = np.asarray(N.T @ a).ravel() / size.replace(0, np.nan).to_numpy()
    if "cyclephase_h" in obs:
        cyc = obs["cyclephase_h"].astype(str).isin(["G1", "S", "G2M"]).to_numpy(float)
        ann["cycling_share"] = np.asarray(N.T @ cyc).ravel() / size.replace(0, np.nan).to_numpy()
    man.write(ann[keep_nh.to_numpy()], "nhood_annotation",
              "Per neighbourhood: size, dominant class and its share, dominant cluster, mean age, cycling share",
              subdir=SUBDIR)
    tot = pd.crosstab(obs[donor_col].astype(str), obs[cls_col].astype(str))
    tot.insert(0, "all_cells", tot.sum(axis=1))
    man.write(tot.reset_index().rename(columns={donor_col: "donor"}), "donor_totals",
              "Cells per donor, overall and per class (abundance denominators)", subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
