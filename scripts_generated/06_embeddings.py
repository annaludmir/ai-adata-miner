#!/usr/bin/env python3
"""06 - Summarise the stored embeddings, and test them for batch mixing.

Both files ship precomputed embeddings (cortex: PCA/UMAP/TSNE/X_scVI,
human_dev: X_Embedding), so the expensive part is already done.  What is
missing is a quantitative read: where each cell class sits, how spread out it
is, how far apart classes are, and -- the one that decides whether the
embedding is trustworthy -- whether neighbourhoods mix donors or segregate them.

The mixing score is a kNN test on a subsample: for each cell, the fraction of
its k nearest neighbours sharing its donor, compared with the fraction expected
if donors were distributed at random.  A ratio near 1 means well mixed; much
greater than 1 means the embedding is organised by batch.

Outputs (csv_exports/<dataset>/06_embeddings/)
  embedding_overview.csv              key, dimensionality, variance per axis
  centroids_<key>_by_<role>.csv       per-group centroid and dispersion
  centroid_distances_<key>_<role>.csv pairwise distances between group centroids
  batch_mixing_<key>.csv              kNN donor/sample mixing scores
  embedding_cells_<key>.csv           subsampled per-cell coords + metadata
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
                          log, read_elem_at, resolve_role)

SCRIPT = "06_embeddings"
SUBDIR = "06_embeddings"
GROUP_ROLES = ["cell_class", "region", "subregion", "cyclephase", "donor", "chemistry"]
MIXING_ROLES = ["donor", "sample", "chemistry"]
SUBSAMPLE_CELLS = 50_000
MIXING_CELLS = 20_000
K_NEIGHBOURS = 30
MAX_LEVELS = 300


def knn_mixing(coords: np.ndarray, labels: pd.Series, k: int) -> pd.DataFrame:
    """Per-label observed vs expected same-label neighbour fraction."""
    from scipy.spatial import cKDTree
    valid = labels.notna().to_numpy()
    coords, labels = coords[valid], labels[valid].astype(str).reset_index(drop=True)
    if len(labels) < k + 1 or labels.nunique() < 2:
        return pd.DataFrame()
    tree = cKDTree(coords)
    # +1 because the first neighbour returned is the point itself
    _, idx = tree.query(coords, k=min(k + 1, len(labels)))
    neigh = idx[:, 1:]
    codes = pd.Categorical(labels).codes
    same = (codes[neigh] == codes[:, None]).mean(axis=1)

    df = pd.DataFrame({"label": labels.to_numpy(), "same_label_fraction": same})
    expected = labels.value_counts(normalize=True)
    out = df.groupby("label").agg(n_cells=("same_label_fraction", "size"),
                                  observed_same_fraction=("same_label_fraction", "mean"))
    out["expected_same_fraction"] = out.index.map(expected).astype(float)
    out["mixing_ratio"] = out["observed_same_fraction"] / out["expected_same_fraction"]
    out["interpretation"] = np.where(
        out["mixing_ratio"] > 3, "strongly segregated (likely batch structure)",
        np.where(out["mixing_ratio"] > 1.5, "partially segregated", "well mixed"))
    return out.reset_index()


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
    keys = list_h5ad_keys(path)
    available = [k for k in config.EMBEDDING_KEYS if k in keys["obsm"]]
    available += [k for k in keys["obsm"] if k not in available]
    if not available:
        log("  no obsm embeddings present")
        man.flush()
        return
    log(f"  obsm keys: {', '.join(keys['obsm'])}")

    roles = {r: resolve_role(obs, r) for r in GROUP_ROLES}
    roles = {r: c for r, c in roles.items()
             if c is not None and obs[c].nunique(dropna=True) <= MAX_LEVELS}
    rng = np.random.default_rng(config.RANDOM_SEED)
    overview = []

    for ekey in available:
        mat = read_elem_at(path, f"obsm/{ekey}")
        if mat is None:
            continue
        if sp.issparse(mat):
            log(f"  skipping obsm/{ekey}: sparse ({mat.shape}) -- this is a "
                "neighbourhood/indicator matrix (e.g. Milo nhoods), not an embedding")
            continue
        mat = np.asarray(mat)
        if mat.ndim != 2 or mat.shape[1] < 2:
            log(f"  skipping obsm/{ekey}: shape {mat.shape} is not a 2-D embedding")
            continue
        if args.limit_cells:
            mat = mat[:args.limit_cells]
        # obsm rows align with the unfiltered cells; apply the same mask.
        if keep is not None and mat.shape[0] == keep.size:
            mat = mat[keep]
        if mat.shape[0] != len(obs):
            log(f"  skipping obsm/{ekey}: {mat.shape[0]} rows != {len(obs)} cells")
            continue
        mat = np.nan_to_num(mat.astype(np.float64), nan=0.0)
        log(f"  {ekey}: {mat.shape[0]:,} x {mat.shape[1]}")
        overview.append({
            "embedding": ekey, "n_cells": mat.shape[0], "n_dims": mat.shape[1],
            **{f"var_dim{i+1}": float(mat[:, i].var()) for i in range(min(5, mat.shape[1]))},
            "total_variance": float(mat.var(axis=0).sum()),
        })

        # Centroids / dispersion use the first few dims so UMAP and a 50-dim PCA
        # stay comparable; distances are reported in that same subspace.
        ndim = min(mat.shape[1], 10)
        sub = mat[:, :ndim]
        for role, col in roles.items():
            lab = obs[col].astype("object")
            valid = lab.notna().to_numpy()
            if valid.sum() == 0:
                continue
            df = pd.DataFrame(sub[valid], columns=[f"dim{i+1}" for i in range(ndim)])
            df["group"] = lab[valid].astype(str).to_numpy()
            cent = df.groupby("group").mean()
            sizes = df.groupby("group").size().rename("n_cells")
            d2 = ((df[cent.columns].to_numpy()
                   - cent.loc[df["group"]].to_numpy()) ** 2).sum(axis=1)
            disp = pd.Series(np.sqrt(d2), index=df.index).groupby(df["group"]).agg(
                dispersion_mean="mean", dispersion_median="median")
            out = cent.join(sizes).join(disp).reset_index()
            out.insert(0, "embedding", ekey)
            man.write(out, f"centroids_{ekey}_by_{role}",
                      f"Centroid and spread of each {col} in obsm/{ekey} (first {ndim} dims)",
                      subdir=SUBDIR)

            if 2 <= len(cent) <= 100:
                c = cent.to_numpy()
                dist = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(axis=2))
                dm = pd.DataFrame(dist, index=cent.index, columns=cent.index)
                dm.index.name = "group"
                man.write(dm.reset_index(), f"centroid_distances_{ekey}_{role}",
                          f"Euclidean distances between {col} centroids in {ekey}",
                          subdir=SUBDIR)

        # -- batch mixing --------------------------------------------------
        n = mat.shape[0]
        take = rng.choice(n, size=min(MIXING_CELLS, n), replace=False)
        take.sort()
        mix_rows = []
        for role in MIXING_ROLES:
            col = resolve_role(obs, role)
            if col is None or obs[col].nunique(dropna=True) < 2:
                continue
            res = knn_mixing(sub[take], obs[col].iloc[take].reset_index(drop=True),
                             K_NEIGHBOURS)
            if res.empty:
                continue
            res.insert(0, "embedding", ekey)
            res.insert(1, "batch_variable", role)
            mix_rows.append(res)
        if mix_rows:
            man.write(pd.concat(mix_rows, ignore_index=True), f"batch_mixing_{ekey}",
                      f"kNN batch-mixing in {ekey}: observed vs expected same-batch "
                      f"neighbour fraction (k={K_NEIGHBOURS}, n={len(take):,} cells)",
                      subdir=SUBDIR)

        # -- subsampled per-cell coordinates, for downstream plotting --------
        take2 = rng.choice(n, size=min(SUBSAMPLE_CELLS, n), replace=False)
        take2.sort()
        cells = pd.DataFrame(mat[take2, :min(mat.shape[1], 3)],
                             columns=[f"{ekey}_{i+1}" for i in range(min(mat.shape[1], 3))])
        cells.insert(0, "cell_id", obs.index[take2].astype(str))
        for role, col in roles.items():
            cells[role] = obs[col].iloc[take2].astype(str).to_numpy()
        if "age_pcw" in obs.columns:
            cells["age_pcw"] = obs["age_pcw"].iloc[take2].to_numpy()
        man.write(cells, f"embedding_cells_{ekey}",
                  f"Subsample of {len(take2):,} cells with {ekey} coordinates and labels",
                  subdir=SUBDIR)

    if overview:
        man.write(pd.DataFrame(overview), "embedding_overview",
                  "Dimensionality and per-axis variance of each stored embedding",
                  subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
