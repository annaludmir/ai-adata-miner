#!/usr/bin/env python3
"""20 - Differentiation pseudotime along radial glia -> IPC -> neuroblast -> neuron.

Developmental age and differentiation state are tangled: an older sample holds
more neurons, and a neuron's expression depends on how far it has
differentiated. A per-cell position along the neurogenic lineage lets step 3
separate the two -- compare young and old donors at the *same* point of
differentiation.

Method (reads obs and one obsm matrix; no expression):
  * cells of the lineage classes (config.LINEAGE_CLASSES), restricted for the
    whole-brain atlas to config.PSEUDOTIME_REGIONS
  * the stored latent space (config.PSEUDOTIME_EMBEDDING: cortex X_scVI,
    human_dev Factors), each dimension standardised over these cells, then
    reduced to config.PSEUDOTIME_PCS principal components
  * a principal path through the class medians in lineage order; each cell is
    projected onto its nearest segment, and pseudotime is the arc length to
    that point, scaled to [0, 1]
  * equal-width bins (config.PSEUDOTIME_BINS) for pseudobulks (script 09)

The path is fixed by the annotated classes, so pseudotime orders cells within
and between classes but cannot discover a branch the labels do not encode.
Run per chemistry like everything else: the path is fitted within a stratum.

Outputs (csv_exports/<dataset>__<chem>/20_pseudotime/)
  cell_pseudotime.csv       obs_row (position in the file), pseudotime, bin
  path_nodes.csv            class medians the path runs through (PC space)
  pseudotime_by_class.csv   pseudotime quantiles per class (ordering check)
  pseudotime_by_class_x_age.csv  median pseudotime per class x age
  pseudotime_bins.csv       cells per bin x class
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import Manifest, load_obs, log, read_elem_at, resolve_role

SCRIPT = "20_lineage_pseudotime"
SUBDIR = "20_pseudotime"


def principal_path(Z: np.ndarray, cls: np.ndarray, order: list[str]) -> tuple[np.ndarray, list[str]]:
    nodes, names = [], []
    for c in order:
        m = cls == c
        if m.sum() >= 10:
            nodes.append(np.median(Z[m], axis=0))
            names.append(c)
    return np.array(nodes), names


def project(Z: np.ndarray, nodes: np.ndarray) -> np.ndarray:
    """Arc-length position of each row of Z on the polyline through nodes, in [0, 1]."""
    seg = np.diff(nodes, axis=0)
    seglen = np.linalg.norm(seg, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    best_d = np.full(len(Z), np.inf)
    best_s = np.zeros(len(Z))
    for k in range(len(seg)):
        v = Z - nodes[k]
        u = np.clip((v @ seg[k]) / max(seglen[k] ** 2, 1e-12), 0.0, 1.0)
        d = np.linalg.norm(v - u[:, None] * seg[k], axis=1)
        better = d < best_d
        best_d[better] = d[better]
        best_s[better] = cum[k] + u[better] * seglen[k]
    return best_s / cum[-1]


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)
    emb_key = config.PSEUDOTIME_EMBEDDING.get(config.base_dataset(key))
    obs, keep = load_obs(path, key, chem, args.limit_cells)
    if obs is None or emb_key is None or len(obs) == 0:
        log("  no chemistry column, no embedding configured, or no cells -- skipping")
        man.flush()
        return
    cls_col = resolve_role(obs, "cell_class")
    if cls_col is None:
        log("  no cell-class column -- skipping")
        man.flush()
        return
    rows = np.nonzero(keep)[0]
    lin = obs[cls_col].astype(str).isin(config.LINEAGE_CLASSES).to_numpy().copy()
    regions = config.PSEUDOTIME_REGIONS.get(config.base_dataset(key))
    reg_col = resolve_role(obs, "region")
    if regions and reg_col is not None:
        lin &= obs[reg_col].astype(str).isin(regions).to_numpy()
    if lin.sum() < 100:
        log(f"  only {int(lin.sum())} lineage cells -- skipping")
        man.flush()
        return
    E = read_elem_at(path, f"obsm/{emb_key}")
    if E is None:
        log(f"  obsm/{emb_key} not present -- skipping")
        man.flush()
        return
    E = np.asarray(E.todense() if hasattr(E, "todense") else E, dtype=np.float64)
    if args.limit_cells:
        E = E[: args.limit_cells]
    E = E[rows[lin]]
    log(f"  {int(lin.sum()):,} lineage cells; obsm/{emb_key} {E.shape[1]} dims"
        + (f"; regions {', '.join(regions)}" if regions and reg_col else ""))
    Z = (E - E.mean(axis=0)) / np.where(E.std(axis=0) > 0, E.std(axis=0), 1.0)
    k = min(config.PSEUDOTIME_PCS, Z.shape[1])
    _, s, Vt = np.linalg.svd(Z - Z.mean(axis=0), full_matrices=False) if len(Z) <= 200_000 else \
        np.linalg.svd(Z[np.random.default_rng(config.RANDOM_SEED).choice(len(Z), 200_000, replace=False)]
                      - Z.mean(axis=0), full_matrices=False)
    P = (Z - Z.mean(axis=0)) @ Vt[:k].T
    cls = obs[cls_col].astype(str).to_numpy()[lin]
    nodes, names = principal_path(P, cls, config.LINEAGE_CLASSES)
    if len(names) < 3:
        log(f"  only {len(names)} lineage classes with cells -- skipping")
        man.flush()
        return
    t = project(P, nodes)
    nb = config.PSEUDOTIME_BINS
    b = np.minimum((t * nb).astype(int), nb - 1)
    cells = pd.DataFrame({"obs_row": rows[lin], "cell_class": cls, "pseudotime": t, "bin": b})
    if "age_pcw" in obs:
        cells["age_pcw"] = obs["age_pcw"].to_numpy()[lin]
    order = {c: i for i, c in enumerate(names)}
    med = cells.groupby("cell_class").pseudotime.median()
    ok = bool(np.all(np.diff([med[c] for c in names]) > 0))
    log(f"  path through {', '.join(names)}; class medians "
        + ", ".join(f"{c} {med[c]:.2f}" for c in names) + (" (ordered)" if ok else " (NOT ordered)"))

    man.write(cells[["obs_row", "pseudotime", "bin"]], "cell_pseudotime",
              f"Per lineage cell: row in the file, pseudotime in [0, 1] and its bin (of {nb})",
              subdir=SUBDIR)
    nodes_df = pd.DataFrame(nodes, columns=[f"PC{i + 1}" for i in range(k)])
    nodes_df.insert(0, "cell_class", names)
    nodes_df["explained_variance_share"] = np.nan
    nodes_df.loc[0, "explained_variance_share"] = float((s[:k] ** 2).sum() / (s ** 2).sum())
    man.write(nodes_df, "path_nodes", f"Class medians defining the path (PCs of standardised {emb_key})",
              subdir=SUBDIR)
    q = cells.groupby("cell_class").pseudotime.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).unstack()
    q.columns = [f"q{int(c * 100)}" for c in q.columns]
    q.insert(0, "n_cells", cells.cell_class.value_counts())
    q["lineage_order"] = [order.get(c, -1) for c in q.index]
    q["path_ordered"] = ok
    man.write(q.sort_values("lineage_order").reset_index(), "pseudotime_by_class",
              "Pseudotime quantiles per class; path_ordered = class medians increase along the lineage",
              subdir=SUBDIR)
    if "age_pcw" in cells:
        ca = cells.groupby(["cell_class", "age_pcw"]).pseudotime.agg(["size", "median", "mean"]).reset_index()
        man.write(ca.rename(columns={"size": "n_cells"}), "pseudotime_by_class_x_age",
                  "Median / mean pseudotime per class x age", subdir=SUBDIR)
    bins = cells.groupby(["bin", "cell_class"]).size().unstack(fill_value=0)
    bins.insert(0, "pseudotime_from", bins.index / nb)
    man.write(bins.reset_index(), "pseudotime_bins", "Cells per pseudotime bin x class", subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
