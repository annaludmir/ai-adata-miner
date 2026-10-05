#!/usr/bin/env python3
"""19 - Radial-glia sub-types: outer vs ventricular radial glia, per cell.

The "Radial glia" label pools ventricular radial glia (vRG, apical, at the
ventricle) and outer radial glia (oRG, basal, in the outer subventricular
zone -- abundant in human cortex from mid-gestation). They divide differently
(oRG undergo mitotic somal translocation), so a change in the radial-glia cell
cycle with age could be a change in sub-type mix rather than in the cells.
This pass scores every radial-glia cell for each programme and exports the
tables needed to tell the two apart.

Per radial-glia cell: mean log1p(CP10K) of oRG markers (HOPX, PTPRZ1, FAM107A,
TNC, MOXD1, LIFR) and of vRG markers (FBXO32, CTGF/CCN2, CYR61/CCN1, PALLD,
PDGFD; old and new symbols both tried). The difference oRG - vRG calls the
cell oRG (> +MARGIN), vRG (< -MARGIN) or ambiguous. Because any margin is a
choice, threshold-free tables by score quintile are exported too.

Streams X once (stage 2). Needs 09's gene selection for the pseudobulk part;
without it, that part is skipped.

Outputs (csv_exports/<dataset>__<chem>/19_rg_subtypes/)
  rg_marker_check.csv                    mean marker expression per called sub-type
  rg_subtype_by_age.csv                  per age: oRG / vRG / ambiguous counts, mean score
  rg_subtype_by_donor.csv                the same per donor
  rg_phase_by_subtype_x_age.csv          phase composition per sub-type x age (if phases)
  rg_phase_by_subtype_x_age_x_depth.csv  ... within UMI quintiles (radial-glia-wide)
  rg_phase_by_score_bin_x_age.csv        phase composition per oRG-vRG score quintile x age
  rg_subtype_x_age__{pseudobulk_counts,group_summary}.csv   pseudobulk per sub-type x age
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

SCRIPT = "19_rg_subtypes"
SUBDIR = "19_rg_subtypes"
RG_LABEL = "Radial glia"
ORG_MARKERS = ["HOPX", "PTPRZ1", "FAM107A", "TNC", "MOXD1", "LIFR"]
VRG_MARKERS = ["FBXO32", "CTGF", "CCN2", "CYR61", "CCN1", "PALLD", "PDGFD"]
MARGIN = 0.2          # log1p(CP10K) units between the two programme scores
N_BINS = 5            # quintiles of the score difference / of UMIs per cell
SUBTYPES = ["oRG", "vRG", "ambiguous"]


def phase_table(df: pd.DataFrame, by: list[str], phase_col: str) -> pd.DataFrame:
    ct = df.groupby(by + [phase_col], observed=True).size().unstack(phase_col, fill_value=0)
    fr = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)
    fr.insert(0, "n_cells", ct.sum(axis=1))
    return fr.reset_index()


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
    class_col = resolve_role(obs, "cell_class")
    if class_col is None:
        log("  no cell-class column; skipping")
        man.flush()
        return
    with XReader(path) as xr:
        n_obs = xr.n_obs
    n_cells = min(args.limit_cells, n_obs) if args.limit_cells else n_obs
    rg = keep & (obs[class_col].astype(str).to_numpy() == RG_LABEL)
    rg[n_cells:] = False
    log(f"  {int(rg.sum()):,} radial-glia cells (chemistry={chem or 'pooled'})")
    if rg.sum() == 0:
        man.flush()
        return

    genes = gene_frame(read_var(path), key)
    sym = genes["symbol"].astype(str).to_numpy()
    o_idx = [i for i, s in enumerate(sym) if s in ORG_MARKERS]
    v_idx = [i for i, s in enumerate(sym) if s in VRG_MARKERS]
    log(f"  oRG markers found: {', '.join(sym[o_idx])} | vRG markers found: {', '.join(sym[v_idx])}")
    if len(o_idx) < 2 or len(v_idx) < 2:
        log("  fewer than two markers for a programme; skipping")
        man.flush()
        return

    # pseudobulk per sub-type x age on 09's exported genes, if available
    sel_path = config.CSV_EXPORTS / ns / "09_pseudobulk" / "gene_selection.csv"
    gene_mask = kept = None
    if sel_path.exists():
        sel = pd.read_csv(sel_path)
        if len(sel) == len(sym):
            gene_mask = sel["selected"].astype(bool).to_numpy()
            kept = sel.loc[gene_mask, "symbol"].astype(str).tolist()
    if gene_mask is None:
        log("  09 gene selection not found -- pseudobulk per sub-type skipped (run 09 first)")
    ages = sorted(obs.loc[rg, "age_pcw"].dropna().unique()) if "age_pcw" in obs else []
    age_pos = {a: i for i, a in enumerate(ages)}
    levels = [f"{s} | {a:g}" for s in SUBTYPES for a in ages]
    agg = (GroupAggregator(len(levels), int(gene_mask.sum()), config.TARGET_SUM)
           if gene_mask is not None and ages else None)

    rows_o, rows_v, rows_u, rows_i = [], [], [], []
    with XReader(path) as xr:
        for start, stop, chunk in xr.iter_chunks(args.chunk_size):
            if start >= n_cells:
                break
            stop = min(stop, n_cells)
            chunk = chunk[: stop - start]
            m = rg[start:stop]
            if not m.any():
                continue
            totals = np.asarray(chunk.sum(axis=1)).ravel()
            sub = chunk[m]
            t = totals[m]
            scale = np.divide(config.TARGET_SUM, t, out=np.zeros_like(t, dtype=float), where=t > 0)
            o = np.log1p(np.asarray(sub[:, o_idx].todense() if hasattr(sub, "todense") else sub[:, o_idx])
                         * scale[:, None]).mean(axis=1)
            v = np.log1p(np.asarray(sub[:, v_idx].todense() if hasattr(sub, "todense") else sub[:, v_idx])
                         * scale[:, None]).mean(axis=1)
            o, v = np.asarray(o).ravel(), np.asarray(v).ravel()
            idx = np.nonzero(m)[0] + start
            rows_o.append(o)
            rows_v.append(v)
            rows_u.append(t)
            rows_i.append(idx)
            if agg is not None:
                diff = o - v
                st = np.where(diff > MARGIN, 0, np.where(diff < -MARGIN, 1, 2))
                a = obs["age_pcw"].to_numpy()[idx]
                codes = np.full(stop - start, -1, dtype=np.int64)
                ok = ~pd.isna(a)
                codes[idx[ok] - start] = st[ok] * len(ages) + np.array([age_pos[x] for x in a[ok]])
                agg.update(chunk, codes, totals_chunk=totals, gene_mask=gene_mask)
    cells = obs.iloc[np.concatenate(rows_i)].copy()
    cells["org_score"] = np.concatenate(rows_o)
    cells["vrg_score"] = np.concatenate(rows_v)
    cells["score_diff"] = cells["org_score"] - cells["vrg_score"]
    cells["total_umis_x"] = np.concatenate(rows_u)
    cells["subtype"] = np.where(cells.score_diff > MARGIN, "oRG",
                                np.where(cells.score_diff < -MARGIN, "vRG", "ambiguous"))
    cells["score_bin"] = pd.qcut(cells.score_diff.rank(method="first"), N_BINS, labels=False)
    cells["depth_bin"] = pd.qcut(cells.total_umis_x.rank(method="first"), N_BINS, labels=False)
    log("  called: " + ", ".join(f"{s} {int((cells.subtype == s).sum()):,}" for s in SUBTYPES))

    # marker sanity: mean expression per called sub-type (from the pass itself)
    chk = cells.groupby("subtype")[["org_score", "vrg_score"]].mean()
    chk.insert(0, "n_cells", cells.subtype.value_counts().reindex(chk.index))
    chk["markers_oRG"] = "|".join(sym[o_idx])
    chk["markers_vRG"] = "|".join(sym[v_idx])
    man.write(chk.reset_index(), "rg_marker_check",
              "Mean oRG and vRG programme scores per called sub-type", subdir=SUBDIR)

    def composition(by: str, label: str) -> None:
        if by is None or by not in cells:
            return
        ct = cells.groupby([by, "subtype"], observed=True).size().unstack("subtype", fill_value=0)
        ct = ct.reindex(columns=SUBTYPES, fill_value=0)
        out = ct.copy()
        out.insert(0, "n_rg", ct.sum(axis=1))
        for s in SUBTYPES:
            out[f"frac_{s}"] = ct[s] / ct.sum(axis=1)
        g = cells.groupby(by, observed=True)
        out["mean_score_diff"] = g["score_diff"].mean()
        out["median_umis"] = g["total_umis_x"].median()
        if label == "donor" and "age_pcw" in cells:
            out["age_pcw"] = g["age_pcw"].first()
        man.write(out.reset_index().rename(columns={by: label}), f"rg_subtype_by_{label}",
                  f"Radial-glia sub-type calls per {label} (margin {MARGIN} log1p CP10K)",
                  subdir=SUBDIR)

    composition("age_pcw", "age")
    composition(resolve_role(obs, "donor"), "donor")

    phase_col = "cyclephase_h" if "cyclephase_h" in cells else None
    if phase_col is not None and "age_pcw" in cells:
        man.write(phase_table(cells, ["subtype", "age_pcw"], phase_col).rename(columns={"age_pcw": "age"}),
                  "rg_phase_by_subtype_x_age", "Phase composition per radial-glia sub-type x age",
                  subdir=SUBDIR)
        man.write(phase_table(cells, ["subtype", "age_pcw", "depth_bin"], phase_col)
                  .rename(columns={"age_pcw": "age"}), "rg_phase_by_subtype_x_age_x_depth",
                  "Phase composition per sub-type x age x UMI quintile (quintiles over all "
                  "radial glia of the stratum)", subdir=SUBDIR)
        man.write(phase_table(cells, ["score_bin", "age_pcw"], phase_col).rename(columns={"age_pcw": "age"}),
                  "rg_phase_by_score_bin_x_age",
                  "Phase composition per quintile of the oRG - vRG score x age (threshold-free; "
                  "bin 4 = most oRG-like)", subdir=SUBDIR)

    if agg is not None:
        man.write(agg.pseudobulk_counts(levels, kept).reset_index(), "rg_subtype_x_age__pseudobulk_counts",
                  "Summed raw counts per gene per radial-glia sub-type x age", subdir=SUBDIR)
        gs = agg.group_summary(levels)
        gs["below_min_cells"] = gs["n_cells"] < config.MIN_CELLS_PER_GROUP
        gs.insert(0, "grouping", "rg_subtype_x_age")
        man.write(gs, "rg_subtype_x_age__group_summary",
                  "Cells and depth per radial-glia sub-type x age", subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
