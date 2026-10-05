#!/usr/bin/env python3
"""09 - Cell-cycle programmes: proliferation over development, and a phase map of genes.

Questions
  A. Within each progenitor type, does proliferation change with age -- the
     cycling fraction, and, among cycling cells, the share in G1 (G1
     lengthening during neurogenesis is a classic expectation), S and G2/M?
  B. Which genes track proliferation across cell clusters, and among cycling
     clusters, which lean to S phase versus G2/M? Summarised per gene list,
     seed NDD panel and co-expression module (08): proliferative or not, and
     which phase.
  C. Do the gene lists that rise with age inside progenitors (06) still rise
     when proliferation-linked genes are set aside? If progenitors cycle less
     later on, genes of non-cycling cells would rise for that reason alone.

Data: cortex carries a per-cell phase call (G1 / S / G2M / Post-M /
Non-cycling); human_dev only a continuous cell-cycle score. So the S-vs-G2/M
parts are cortex only, and human_dev contributes proliferation association.
The phase calls come from cell-cycle marker expression, so the canonical
markers' own phase assignment is partly circular -- it serves as a sanity
check, and the informative part is everything else.

Method
  A. Per stratum x cell class, the per-age value (cycling fraction, mean cycle
     score, and G1/S/G2M shares among cycling cells where phase fractions per
     class x age exist) vs age: Spearman, exact permutation, v2/v3 signed
     Stouffer, tiered as elsewhere. The phase shares are also computed at
     matched depth -- within UMI quintiles, then averaged with fixed weights --
     because phase calls come from marker expression and UMIs per cell fall
     with age in cortex; a trend that vanishes there is a depth artefact.
  B. Cluster pseudobulks as in 07/08. Proliferation association = Spearman of
     each gene's log2 TMM-CPM with the cluster's cycling fraction (cortex) or
     mean cycle score (human_dev). Phase bias (cortex) = partial Spearman with
     the cluster's S share of S+G2M cells, controlling for its cycling
     fraction, over clusters with >= 5% of cells in S or G2/M. Per set: mean
     association vs random sets matched on level x spread, v2/v3 combined.
  C. 06's list-level age test (mean rho with age vs matched random genes)
     repeated on the list's genes with |proliferation association| < 0.2 in
     both chemistries.

Inputs (csv_exports/):
  <ds>__<chem>/03_cellcycle/{proliferation_trajectory,phase_fractions_by_cell_class_x_age,
                             phase_fractions_by_cluster_<c>,cycle_scores_by_cluster_<c>}.csv
  <ds>__<chem>/09_pseudobulk/<finest clustering>__pseudobulk_counts.csv (+ 04_clusters profile)
  <ds>__<chem>/11_panels/panel_coverage.csv
  results/03_age_trends_within_cell_class/age_trends_per_stratum.csv
  results/08_coexpression_modules/modules.csv (optional)
  gene lists: config.GENE_LISTS_DIR (optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "09_cell_cycle_programs"
TITLE = "Cell-cycle programmes: proliferation over development and a phase map of genes"
MIN_AGES = 5
MIN_POINT_CELLS = 50
MIN_CYCLING_CELLS = 20       # cycling cells needed at an age point to compute phase shares
MIN_BIN_CYCLING = 5          # ... and per depth quintile for depth-matched shares
MIN_DEPTH_WEIGHT = 0.6       # depth quintiles present must carry >= this share of the weight
MIN_CYCLING_SHARE = 0.05     # clusters with >= this fraction of cells in S or G2/M
PROLIF_R = 0.3               # |rho| for a gene to count as proliferation-linked
PHASE_R = 0.2                # |partial rho| for an S or G2/M lean
INDEPENDENT_R = 0.2          # |rho| below this in both chemistries = cycle-independent
N_RANDOM = 2000
N_BINS = 5
PROGENITORS = ("Radial glia", "Neuronal IPC", "Glioblast")
MAX_LISTED = 8


# ---------------------------------------------------------------------------
# A. proliferation over development
# ---------------------------------------------------------------------------
def depth_standardised_shares(dp: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """G1/S/G2M shares of cycling cells per class x age, compared at matched depth.

    Phase calls come from marker expression and can drift with UMIs per cell,
    which falls with age in cortex. Within each UMI quintile (quintiles are
    per cell class, over all ages) the shares are computed separately, then
    averaged with fixed weights -- each quintile's share of the class's cycling
    cells over all ages -- so every age point is read at the same depth mix.
    An age point needs quintiles covering >= MIN_DEPTH_WEIGHT of that weight,
    each with >= MIN_BIN_CYCLING cycling cells. Also returns median UMIs per
    cell per class x age, to show the depth drift itself.
    """
    if not {"G1", "S", "G2M"} <= set(dp.columns):
        return {}
    dp = dp.copy()
    for ph in ("G1", "S", "G2M"):
        dp[f"n_{ph}"] = dp[ph].fillna(0) * dp.n_cells
    dp["n_cyc"] = dp[["n_G1", "n_S", "n_G2M"]].sum(axis=1)
    out = {f"{ph} share of cycling cells (depth-matched)": [] for ph in ("G1", "S", "G2M")}
    out["median UMIs per cell"] = []
    for cls, g in dp.groupby("cell_class"):
        w = g.groupby("depth_bin")["n_cyc"].sum()
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        for age, h in g.groupby("age"):
            out["median UMIs per cell"].append(
                {"cell_class": cls, "age_pcw": age,
                 "value": float(np.average(h.median_umis, weights=h.n_cells))})
            h = h[h.n_cyc >= MIN_BIN_CYCLING].set_index("depth_bin")
            cover = w.reindex(h.index).fillna(0)
            if cover.sum() < MIN_DEPTH_WEIGHT:
                continue
            for ph in ("G1", "S", "G2M"):
                share = h[f"n_{ph}"] / h["n_cyc"]
                out[f"{ph} share of cycling cells (depth-matched)"].append(
                    {"cell_class": cls, "age_pcw": age,
                     "value": float((share * cover).sum() / cover.sum())})
    return {k: pd.DataFrame(v, columns=["cell_class", "age_pcw", "value"]) for k, v in out.items() if v}


def trajectories(out: C.Output) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        out.used(f"{n}/03_cellcycle/proliferation_trajectory.csv")
        tr = C.csv(n, "03_cellcycle/proliferation_trajectory.csv")
        tr = tr[tr.n_cells >= MIN_POINT_CELLS]
        metrics = {}
        for col, label in (("fraction_cycling", "fraction cycling"),
                           ("cell_cycle_score_mean", "mean cycle score")):
            if col in tr and tr[col].notna().any():
                metrics[label] = tr[["cell_class", "age_pcw", col]].rename(columns={col: "value"})
        pf_path = C.EXPORTS / n / "03_cellcycle" / "phase_fractions_by_cell_class_x_age.csv"
        if pf_path.exists():
            out.used(f"{n}/03_cellcycle/phase_fractions_by_cell_class_x_age.csv")
            pf = pd.read_csv(pf_path)
            pf = pf[pf.n_cells >= MIN_POINT_CELLS]
            if {"G1", "S", "G2M"} <= set(pf.columns):
                cyc = pf[["G1", "S", "G2M"]].sum(axis=1)
                ok = cyc * pf.n_cells >= MIN_CYCLING_CELLS   # enough cycling cells for shares
                for ph in ("G1", "S", "G2M"):
                    d = pf.loc[ok, ["cell_class", "age"]].rename(columns={"age": "age_pcw"})
                    d["value"] = (pf.loc[ok, ph] / cyc[ok]).to_numpy()
                    metrics[f"{ph} share of cycling cells"] = d
        dp_path = C.EXPORTS / n / "03_cellcycle" / "phase_fractions_by_cell_class_x_age_x_depth.csv"
        if dp_path.exists():
            out.used(f"{n}/03_cellcycle/phase_fractions_by_cell_class_x_age_x_depth.csv")
            metrics.update(depth_standardised_shares(pd.read_csv(dp_path)))
        for label, d in metrics.items():
            for cls, g in d.dropna().groupby("cell_class"):
                g = g.sort_values("age_pcw")
                if len(g) < MIN_AGES or g["value"].std() == 0:
                    continue
                ages, vals = g.age_pcw.to_numpy(float), g.value.to_numpy(float)
                rho = C.spearman_rows(vals[None, :], ages)
                p, exact = C.spearman_perm_p(rho, ages, X=vals[None, :])
                rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "metric": label,
                             "n_ages": len(g), "first": vals[0], "last": vals[-1],
                             "age_range": f"{ages.min():g}-{ages.max():g}",
                             "rho": float(rho[0]), "perm_p": float(p[0]), "exact": exact})
    per = pd.DataFrame(rows)
    comb = []
    for (ds, cls, metric), g in per.groupby(["dataset", "cell_class", "metric"], sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        Z, pc = C.signed_stouffer([np.array([a.rho]), np.array([b.rho])],
                                  [np.array([a.perm_p]), np.array([b.perm_p])],
                                  [np.sqrt(a.n_ages), np.sqrt(b.n_ages)])
        comb.append({"dataset": ds, "cell_class": cls, "metric": metric,
                     "rho_v2": a.rho, "p_v2": a.perm_p, "range_v2": f"{a['first']:.3g}->{a['last']:.3g}",
                     "rho_v3": b.rho, "p_v3": b.perm_p, "range_v3": f"{b['first']:.3g}->{b['last']:.3g}",
                     "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
    comb = pd.DataFrame(comb)
    if not comb.empty:
        comb["combined_q"] = np.nan
        for _, ix in comb.groupby("dataset").groups.items():
            comb.loc[ix, "combined_q"] = C.bh(comb.loc[ix, "combined_p"])
        comb["tier"] = C.replication_tier(comb.rho_v2, comb.p_v2, comb.rho_v3, comb.p_v3,
                                          comb.combined_q)
        comb["direction"] = np.where(comb.stouffer_z > 0, "rises", "falls")
        comb = comb.sort_values(["dataset", "combined_p"])
    out.write(per, "proliferation_trajectories_per_stratum",
              "Per stratum x class: proliferation metric vs age (exact permutation p)")
    out.write(comb, "proliferation_trajectories_combined",
              "v2 x v3 combined proliferation trends per class and metric; tier replicated / supported")
    return per, comb


# ---------------------------------------------------------------------------
# B. gene phase map
# ---------------------------------------------------------------------------
def cluster_covariates(n: str, cols: list[str]) -> pd.DataFrame:
    """Per cluster: cycling fraction, S share of S+G2M (cortex), mean cycle score."""
    ds = n.split("__")[0]
    clustering = C.COEXPR_GROUPING[ds].replace("cluster_", "", 1)
    cov = pd.DataFrame(index=cols)
    pf_path = C.EXPORTS / n / "03_cellcycle" / f"phase_fractions_by_cluster_{clustering}.csv"
    if pf_path.exists():
        pf = pd.read_csv(pf_path)
        pf = pf.set_index(pf.columns[0])
        pf.index = pf.index.astype(str)
        pf = pf.reindex(cols)
        if {"G1", "S", "G2M"} <= set(pf.columns):
            cov["cycling"] = pf[["G1", "S", "G2M"]].sum(axis=1)
            sg = pf["S"] + pf["G2M"]
            cov["s_share"] = np.where(sg >= MIN_CYCLING_SHARE, pf["S"] / sg.replace(0, np.nan), np.nan)
    cs_path = C.EXPORTS / n / "03_cellcycle" / f"cycle_scores_by_cluster_{clustering}.csv"
    if cs_path.exists():
        cs = pd.read_csv(cs_path)
        cs["group_level"] = cs["group_level"].astype(str)
        cov["cycle_score"] = cs.set_index("group_level")["cell_cycle_score_mean"].reindex(cols)
    return cov


def partial_spearman(X: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Spearman of each row of X with y, controlling for z (ranks residualised on z)."""
    R = C.rank_rows(X)
    ry, rz = C.stats.rankdata(y), C.stats.rankdata(z)
    R = R - R.mean(axis=1, keepdims=True)
    ry, rz = ry - ry.mean(), rz - rz.mean()
    zz = rz @ rz
    Rres = R - np.outer(R @ rz / zz, rz)
    yres = ry - (ry @ rz / zz) * rz
    den = np.linalg.norm(Rres, axis=1) * np.linalg.norm(yres)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, Rres @ yres / den, np.nan)


def phase_map(out: C.Output) -> tuple[pd.DataFrame, dict]:
    per, mats = [], {}
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        clustering = C.COEXPR_GROUPING[ds].replace("cluster_", "", 1)
        lc, _ = C.cluster_expression(n)
        cov = cluster_covariates(n, list(lc.columns))
        for f in (f"phase_fractions_by_cluster_{clustering}", f"cycle_scores_by_cluster_{clustering}"):
            if (C.EXPORTS / n / "03_cellcycle" / f"{f}.csv").exists():
                out.used(f"{n}/03_cellcycle/{f}.csv")
        prolif_col = "cycling" if "cycling" in cov and cov["cycling"].notna().sum() > 20 else "cycle_score"
        ok = cov[prolif_col].notna().to_numpy()
        X = lc.to_numpy(float)
        rho_p = C.spearman_rows(X[:, ok], cov[prolif_col].to_numpy(float)[ok])
        d = pd.DataFrame({"dataset": ds, "chemistry": chem, "gene": lc.index,
                          "proliferation_measure": prolif_col, "rho_proliferation": rho_p,
                          "mean_log2cpm": X.mean(axis=1)})
        if "s_share" in cov and cov["s_share"].notna().sum() >= 20:
            m = cov["s_share"].notna().to_numpy()
            d["rho_s_vs_g2m"] = partial_spearman(X[:, m], cov["s_share"].to_numpy(float)[m],
                                                 cov["cycling"].to_numpy(float)[m])
            d["n_cycling_clusters"] = int(m.sum())
        per.append(d)
        mats[(ds, chem)] = lc
    per = pd.concat(per, ignore_index=True)
    out.write(per, "gene_phase_map_per_stratum",
              "Per gene x stratum: Spearman with cluster proliferation; partial Spearman with the "
              "S share of S+G2M cells (cortex), controlling for proliferation")
    # combine chemistries per dataset and classify
    rows = []
    for ds, g in per.groupby("dataset"):
        w = g.pivot_table(index="gene", columns="chemistry",
                          values=[c for c in ("rho_proliferation", "rho_s_vs_g2m") if c in g])
        w.columns = [f"{a}_{b}" for a, b in w.columns]
        w = w.dropna(subset=["rho_proliferation_v2", "rho_proliferation_v3"])
        p2, p3 = w.rho_proliferation_v2, w.rho_proliferation_v3
        cls = np.where((p2 >= PROLIF_R) & (p3 >= PROLIF_R), "proliferative",
                       np.where((p2 <= -PROLIF_R) & (p3 <= -PROLIF_R), "anti-proliferative",
                                np.where((p2.abs() < INDEPENDENT_R) & (p3.abs() < INDEPENDENT_R),
                                         "cycle-independent", "unresolved")))
        w["proliferation_class"] = cls
        if {"rho_s_vs_g2m_v2", "rho_s_vs_g2m_v3"} <= set(w.columns):
            s2, s3 = w.rho_s_vs_g2m_v2, w.rho_s_vs_g2m_v3
            lean = np.where((s2 >= PHASE_R) & (s3 >= PHASE_R), "S",
                            np.where((s2 <= -PHASE_R) & (s3 <= -PHASE_R), "G2/M", ""))
            w["phase_lean"] = np.where(w.proliferation_class == "proliferative", lean, "")
        else:
            w["phase_lean"] = ""
        w.insert(0, "dataset", ds)
        rows.append(w.reset_index())
    comb = pd.concat(rows, ignore_index=True)
    out.write(comb, "gene_phase_map",
              f"Per gene: proliferation class (|rho| >= {PROLIF_R} in both chemistries) and, for "
              f"proliferative genes, S or G2/M lean (|partial rho| >= {PHASE_R} in both; cortex)")
    return comb, mats


def set_phase_profile(out: C.Output, gene_map: pd.DataFrame, mats: dict,
                      sets_by_ds: dict, rng: np.random.Generator) -> pd.DataFrame:
    """Per gene set: mean proliferation association and S-vs-G2/M lean vs matched null."""
    per_stratum = pd.read_csv(C.RESULTS / SLUG / "gene_phase_map_per_stratum.csv")
    rows = []
    for (ds, chem), lc in mats.items():
        g = per_stratum[(per_stratum.dataset == ds) & (per_stratum.chemistry == chem)].set_index("gene")
        g = g.loc[lc.index]
        bins = C.expression_bins(lc, N_BINS)
        pools = {b: np.nonzero(bins == b)[0] for b in np.unique(bins)}
        pos = {x: i for i, x in enumerate(lc.index)}
        for sname, genes in sets_by_ds[ds].items():
            idx = np.array([pos[x] for x in genes if x in pos])
            if idx.size < 5:
                continue
            rand = np.stack([rng.choice(pools[bins[i]], size=N_RANDOM) for i in idx], axis=1)
            rec = {"dataset": ds, "chemistry": chem, "gene_set": sname, "n_genes": int(idx.size)}
            for col in ("rho_proliferation", "rho_s_vs_g2m"):
                if col not in g or g[col].isna().all():
                    continue
                v = g[col].to_numpy(float)
                ok_idx = idx[np.isfinite(v[idx])]
                if ok_idx.size < 5:
                    continue
                obs = float(np.nanmean(v[ok_idx]))
                vv = np.nan_to_num(v, nan=0.0)
                null = C.set_mean_rows(vv[:, None], rand)[:, 0]
                sd = null.std()
                rec[f"{col}_mean"] = obs
                rec[f"{col}_null"] = float(null.mean())
                rec[f"{col}_effect"] = (obs - null.mean()) / sd if sd > 0 else np.nan
                rec[f"{col}_p"] = float((np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1)
                                        / (N_RANDOM + 1))
            rows.append(rec)
    per = pd.DataFrame(rows)
    comb = []
    for (ds, sname), g in per.groupby(["dataset", "gene_set"], sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        rec = {"dataset": ds, "gene_set": sname,
               "n_genes": int(g.n_genes.min())}
        for col in ("rho_proliferation", "rho_s_vs_g2m"):
            if f"{col}_effect" not in g or g[f"{col}_effect"].isna().any():
                continue
            a, b = g.loc["v2"], g.loc["v3"]
            Z, pc = C.signed_stouffer([np.array([a[f"{col}_effect"]]), np.array([b[f"{col}_effect"]])],
                                      [np.array([a[f"{col}_p"]]), np.array([b[f"{col}_p"]])], [1.0, 1.0])
            rec.update({f"{col}_v2": a[f"{col}_mean"], f"{col}_v3": b[f"{col}_mean"],
                        f"{col}_null_v2": a[f"{col}_null"], f"{col}_null_v3": b[f"{col}_null"],
                        f"{col}_effect_v2": a[f"{col}_effect"], f"{col}_effect_v3": b[f"{col}_effect"],
                        f"{col}_p_v2": a[f"{col}_p"], f"{col}_p_v3": b[f"{col}_p"],
                        f"{col}_combined_p": float(pc[0])})
        gm = gene_map[(gene_map.dataset == ds) & gene_map.gene.isin(sets_by_ds[ds].get(sname, []))]
        rec["n_proliferative"] = int((gm.proliferation_class == "proliferative").sum())
        rec["n_S_lean"] = int((gm.phase_lean == "S").sum())
        rec["n_G2M_lean"] = int((gm.phase_lean == "G2/M").sum())
        rec["n_anti_proliferative"] = int((gm.proliferation_class == "anti-proliferative").sum())
        rec["n_cycle_independent"] = int((gm.proliferation_class == "cycle-independent").sum())
        rec["S_genes"] = "|".join(gm.loc[gm.phase_lean == "S", "gene"])
        rec["G2M_genes"] = "|".join(gm.loc[gm.phase_lean == "G2/M", "gene"])
        comb.append(rec)
    comb = pd.DataFrame(comb)
    # (the S-vs-G2/M columns are only interpretable for sets that are proliferative)
    for col in ("rho_proliferation", "rho_s_vs_g2m"):
        if f"{col}_combined_p" in comb:
            comb[f"{col}_q"] = np.nan
            for _, ix in comb.groupby("dataset").groups.items():
                ok = comb.loc[ix, f"{col}_combined_p"].notna()
                comb.loc[ok[ok].index, f"{col}_q"] = C.bh(comb.loc[ok[ok].index, f"{col}_combined_p"])
            comb[f"{col}_tier"] = C.replication_tier(
                comb[f"{col}_effect_v2"], comb[f"{col}_p_v2"], comb[f"{col}_effect_v3"],
                comb[f"{col}_p_v3"], comb[f"{col}_q"].fillna(1.0))
    out.write(per, "set_phase_profile_per_stratum",
              "Per gene set x stratum: mean proliferation association / S-vs-G2M lean vs matched null")
    out.write(comb, "set_phase_profile",
              "Per gene set: proliferation association and S-vs-G2/M lean (v2/v3 combined, tiered), "
              "with counts of proliferative, S-, G2/M-leaning, anti-proliferative genes")
    return comb


# ---------------------------------------------------------------------------
# C. age trends in progenitors, with proliferation-linked genes set aside
# ---------------------------------------------------------------------------
def age_trends_without_proliferation(out: C.Output, gene_map: pd.DataFrame, lists: dict,
                                     rng: np.random.Generator) -> pd.DataFrame:
    path = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_per_stratum.csv"
    per = pd.read_csv(C.require(path), low_memory=False, dtype={"panels": str})
    out.used("results/03_age_trends_within_cell_class/age_trends_per_stratum.csv")
    rows = []
    for (ds, chem, cls), g in per.groupby(["dataset", "chemistry", "cell_class"], sort=False):
        if cls not in PROGENITORS:
            continue
        gm = gene_map[gene_map.dataset == ds].set_index("gene")
        g = g.reset_index(drop=True)
        indep = g.gene.map(gm["proliferation_class"]).eq("cycle-independent").to_numpy()
        prol = g.gene.map(gm["rho_proliferation_v2" if chem == "v2" else "rho_proliferation_v3"])
        rho = g.spearman_rho.to_numpy(float)
        both = np.isfinite(prol.to_numpy(float))
        r_all = float(pd.Series(rho[both]).corr(pd.Series(prol.to_numpy(float)[both]), method="spearman"))
        level = g.mean_log2cpm.to_numpy(float)
        bins = np.digitize(level, np.quantile(level, np.linspace(0, 1, 11)[1:-1]))
        pools = [np.nonzero((bins == b) & indep)[0] for b in range(10)]
        pos = {x: i for i, x in enumerate(g.gene)}
        for name, genes in lists.get(ds, {}).items():
            idx = np.array([pos[x] for x in genes if x in pos and indep[pos[x]]])
            if idx.size < 5 or any(pools[bins[i]].size == 0 for i in idx):
                continue
            obs = float(rho[idx].mean())
            rand = np.stack([rng.choice(pools[bins[i]], size=N_RANDOM) for i in idx], axis=1)
            null = C.set_mean_rows(rho[:, None], rand)[:, 0]
            sd = null.std()
            rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "gene_list": name,
                         "n_cycle_independent_genes": int(idx.size),
                         "n_ages": int(g.n_ages.iloc[0]),
                         "mean_rho_age": obs, "null_mean": float(null.mean()),
                         "effect_vs_null_sd": (obs - null.mean()) / sd if sd > 0 else np.nan,
                         "perm_p": float((np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1)
                                         / (N_RANDOM + 1)),
                         "class_wide_rho_age_vs_proliferation": r_all})
    perstr = pd.DataFrame(rows)
    comb = []
    for (ds, cls, name), g in perstr.groupby(["dataset", "cell_class", "gene_list"], sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        Z, pc = C.signed_stouffer([np.array([a.effect_vs_null_sd]), np.array([b.effect_vs_null_sd])],
                                  [np.array([a.perm_p]), np.array([b.perm_p])],
                                  [np.sqrt(a.n_ages), np.sqrt(b.n_ages)])
        comb.append({"dataset": ds, "cell_class": cls, "gene_list": name,
                     "n_genes_v2": int(a.n_cycle_independent_genes), "n_genes_v3": int(b.n_cycle_independent_genes),
                     "mean_rho_v2": a.mean_rho_age, "effect_v2": a.effect_vs_null_sd, "p_v2": a.perm_p,
                     "mean_rho_v3": b.mean_rho_age, "effect_v3": b.effect_vs_null_sd, "p_v3": b.perm_p,
                     "class_wide_rho_age_vs_prolif_v2": a.class_wide_rho_age_vs_proliferation,
                     "class_wide_rho_age_vs_prolif_v3": b.class_wide_rho_age_vs_proliferation,
                     "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
    comb = pd.DataFrame(comb)
    if not comb.empty:
        comb["combined_q"] = np.nan
        for _, ix in comb.groupby("dataset").groups.items():
            comb.loc[ix, "combined_q"] = C.bh(comb.loc[ix, "combined_p"])
        comb["tier"] = C.replication_tier(comb.effect_v2, comb.p_v2, comb.effect_v3, comb.p_v3,
                                          comb.combined_q)
        comb["direction"] = np.where(comb.stouffer_z > 0, "rises", "falls")
        comb = comb.sort_values(["dataset", "combined_p"])
    out.write(perstr, "age_trends_cycle_independent_per_stratum",
              "List age trend in progenitors using only cycle-independent genes, per stratum")
    out.write(comb, "age_trends_cycle_independent",
              "v2 x v3 combined list age trends in progenitors, cycle-independent genes only")
    return comb


# ---------------------------------------------------------------------------
def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)

    _, traj = trajectories(out)
    gene_map, mats = phase_map(out)

    # gene sets: user lists, seed NDD + cell-cycle panels, robust 08 modules
    lists = {ds: C.mapped_lists(ds) for ds in C.DATASETS} if C.gene_lists() else {}
    sets_by_ds = {}
    for ds in C.DATASETS:
        s = {f"list:{k}": v for k, v in lists.get(ds, {}).items()}
        n = C.ns(ds, "v2")
        out.used(f"{n}/11_panels/panel_coverage.csv")
        pan = C.panels(n)
        for (grp, pname), g in pan[pan.panel_group.isin(["ndd", "cell_cycle"])].groupby(["panel_group", "panel"]):
            s[f"seed:{pname}"] = list(g.gene)
        sets_by_ds[ds] = s
    mod_path = C.RESULTS / "08_coexpression_modules" / "modules.csv"
    if mod_path.exists():
        out.used("results/08_coexpression_modules/modules.csv")
        mods = pd.read_csv(mod_path)
        for r in mods[mods.robust].itertuples():
            sets_by_ds[r.dataset][f"module:{r.module}"] = r.genes.split("|")
    prof = set_phase_profile(out, gene_map, mats, sets_by_ds, rng)
    cind = age_trends_without_proliferation(out, gene_map, lists, rng) if lists else pd.DataFrame()

    figures(out, traj, prof)

    # ---- findings -------------------------------------------------------------
    f = []
    if not traj.empty:
        for ds in C.DATASETS:
            t = traj[(traj.dataset == ds) & traj.cell_class.isin(PROGENITORS) & (traj.tier != "")]
            if t.empty:
                f.append(f"**{ds}: no replicated proliferation trend in progenitors.**")
                continue
            f.append(f"**{ds}: proliferation over development in progenitors** -- "
                     + "; ".join(f"{r.metric} {r.direction} with age in {r.cell_class} "
                                 f"({r.range_v2} in v2, {r.range_v3} in v3; {r.tier})"
                                 for r in t.itertuples()) + ".")
        g1 = traj[traj.metric.str.startswith("G1 share") & traj.cell_class.isin(PROGENITORS)]
        if len(g1):
            rep = g1[g1.tier == "replicated"]
            f.append("**G1 lengthening** (G1 share among cycling progenitors rising with age): "
                     + ("; ".join(f"{r.dataset} {r.cell_class} {r.direction} ({r.range_v2} / {r.range_v3})"
                                  for r in rep.itertuples()) if len(rep) else
                        "not replicated in either chemistry pair -- "
                        + "; ".join(f"{r.dataset} {r.cell_class}: rho {r.rho_v2:+.2f} / {r.rho_v3:+.2f}"
                                    for r in g1.itertuples()))
                     + ".")
        dm = traj[traj.metric.str.contains("depth-matched") & traj.cell_class.isin(PROGENITORS)]
        if len(dm):
            raw = traj.set_index(["dataset", "cell_class", "metric"])
            parts = []
            for r in dm.itertuples():
                base = r.metric.replace(" (depth-matched)", "")
                rr = raw.loc[(r.dataset, r.cell_class, base)] if (r.dataset, r.cell_class, base) in raw.index else None
                if rr is None or (rr.tier == "" and r.tier == ""):
                    continue
                parts.append(f"{r.dataset} {r.cell_class} {base.split(' share')[0]} share: raw rho "
                             f"{rr.rho_v2:+.2f}/{rr.rho_v3:+.2f} ({rr.tier or 'n.s.'}), depth-matched "
                             f"{r.rho_v2:+.2f}/{r.rho_v3:+.2f} ({r.tier or 'n.s.'})")
            umi = traj[(traj.metric == "median UMIs per cell") & traj.cell_class.isin(PROGENITORS)]
            f.append("**Depth check on phase shares** (cells compared within UMI quintiles; "
                     "trends that vanish there are likely depth artefacts of the phase calls): "
                     + ("; ".join(parts) if parts else "no phase-share trend to check")
                     + ". UMIs per cell vs age: "
                     + ", ".join(f"{r.dataset} {r.cell_class} rho {r.rho_v2:+.2f}/{r.rho_v3:+.2f}"
                                 for r in umi.itertuples()) + ".")
        if not any(m.startswith("G1 share") for m in traj.metric.unique()):
            f.append("**G1 lengthening not testable yet**: phase fractions per cell class x age are "
                     "exported by script 03 from this version on; re-run stage 1 to add them.")
    for ds in C.DATASETS:
        gm = gene_map[gene_map.dataset == ds]
        cnt = gm.proliferation_class.value_counts()
        lean = gm.phase_lean.value_counts()
        f.append(f"**{ds} gene phase map** ({len(gm):,} genes): "
                 f"{cnt.get('proliferative', 0)} proliferative"
                 + (f" ({lean.get('S', 0)} lean S, {lean.get('G2/M', 0)} lean G2/M)" if lean.get('S', 0) + lean.get('G2/M', 0) else "")
                 + f", {cnt.get('anti-proliferative', 0)} anti-proliferative, "
                 f"{cnt.get('cycle-independent', 0)} cycle-independent.")
    if not prof.empty:
        chk = prof[prof.gene_set.isin(["seed:s_phase", "seed:g2m_phase"])]
        if len(chk) and "rho_s_vs_g2m_effect_v2" in chk:
            f.append("**Sanity check against the seed phase panels** (partly circular: phase calls "
                     "come from such genes): "
                     + "; ".join(f"{r.dataset} {r.gene_set.split(':')[1]}: lean "
                                 f"{r.rho_s_vs_g2m_effect_v2:+.1f}/{r.rho_s_vs_g2m_effect_v3:+.1f} null SDs "
                                 f"(+ = S), {r.n_S_lean} S / {r.n_G2M_lean} G2/M genes"
                                 for r in chk.dropna(subset=["rho_s_vs_g2m_effect_v2"]).itertuples()) + ".")
        for ds in C.DATASETS:
            p = prof[(prof.dataset == ds) & ~prof.gene_set.str.startswith("module:")]
            if p.empty:
                continue
            parts = []
            for r in p.itertuples():
                if getattr(r, "rho_proliferation_tier", "") != "replicated":
                    continue
                sign = "proliferative" if r.rho_proliferation_effect_v2 > 0 else "anti-proliferative"
                lean = ""
                # a phase lean only means something for genes expressed in cycling cells
                if sign == "proliferative" and getattr(r, "rho_s_vs_g2m_tier", "") == "replicated":
                    lean = ", leans S" if r.rho_s_vs_g2m_effect_v2 > 0 else ", leans G2/M"
                comp = (f"; {r.n_S_lean} S-lean, {r.n_G2M_lean} G2/M-lean genes"
                        if r.n_S_lean + r.n_G2M_lean else "")
                parts.append(f"{r.gene_set} {sign}{lean} (rho {r.rho_proliferation_v2:+.2f}/"
                             f"{r.rho_proliferation_v3:+.2f} vs {r.rho_proliferation_null_v2:+.2f} "
                             f"matched{comp})")
            f.append(f"**{ds}: gene sets tied to proliferation, replicated** -- "
                     + ("; ".join(parts) if parts else "none") + ".")
            mp = prof[(prof.dataset == ds) & prof.gene_set.str.startswith("module:")]
            if "rho_s_vs_g2m_tier" in mp and "rho_proliferation_tier" in mp:
                lean = mp[(mp.rho_s_vs_g2m_tier == "replicated") & (mp.rho_proliferation_tier == "replicated")
                          & (mp.rho_proliferation_effect_v2 > 0)]
                if len(lean):
                    f.append(f"**{ds}: proliferative modules with a replicated phase lean** -- "
                             + "; ".join(f"{r.gene_set.split(':')[1]} leans "
                                         f"{'S' if r.rho_s_vs_g2m_effect_v2 > 0 else 'G2/M'}"
                                         for r in lean.itertuples()) + ".")
    if not cind.empty:
        for ds in C.DATASETS:
            c = cind[cind.dataset == ds]
            if c.empty:
                continue
            rep = c[c.tier == "replicated"]
            wide = c.drop_duplicates("cell_class")
            f.append(f"**{ds}: list age trends in progenitors with proliferation-linked genes set aside** "
                     + ("(region-confounded, see 03) " if ds == "human_dev" else "") + "-- "
                     + ("; ".join(f"{r.gene_list} still {r.direction} in {r.cell_class} "
                                  f"({r.n_genes_v2}/{r.n_genes_v3} cycle-independent genes, mean rho "
                                  f"{r.mean_rho_v2:+.2f}/{r.mean_rho_v3:+.2f})"
                                  for r in rep.head(MAX_LISTED).itertuples())
                        if len(rep) else "no list trend survives")
                     + (f" (+{len(rep) - MAX_LISTED} more)" if len(rep) > MAX_LISTED else "")
                     + ". Class-wide, genes' age trends vs their proliferation association: "
                     + ", ".join(f"{r.cell_class} rho {r.class_wide_rho_age_vs_prolif_v2:+.2f}/"
                                 f"{r.class_wide_rho_age_vs_prolif_v3:+.2f}" for r in wide.itertuples())
                     + " (negative = genes of cycling cells fall with age).")

    out.summary(
        TITLE,
        "How does proliferation change over development within progenitor types; which genes "
        "and gene groups follow proliferation, and which lean to S or G2/M; and do list-level age "
        "trends in progenitors survive once proliferation-linked genes are set aside?",
        ["A: per stratum x cell class, cycling fraction (cortex; cells not in Non-cycling or "
         "Post-M), mean cell-cycle score, and G1/S/G2M shares among cycling cells (where 03 "
         f"exports phase fractions per class x age), vs age over points with >= {MIN_POINT_CELLS} "
         "cells; Spearman, exact permutation p, v2/v3 signed Stouffer, BH per dataset, tiered.",
         "B: cluster pseudobulks as in 07/08 (log2 TMM-CPM). Proliferation association = "
         "Spearman with the cluster's cycling fraction (cortex) or mean cycle score (human_dev). "
         "S-vs-G2/M lean (cortex) = partial Spearman with S/(S+G2M) controlling for cycling "
         f"fraction, over clusters with >= {MIN_CYCLING_SHARE:.0%} of cells in S or G2/M. Gene "
         f"classes need |rho| >= {PROLIF_R} (lean: >= {PHASE_R}) in both chemistries; "
         f"cycle-independent = |rho| < {INDEPENDENT_R} in both.",
         f"Set profiles: mean association of a set's genes vs {N_RANDOM:,} random sets matched "
         "member by member on mean level x spread; v2/v3 combined, BH per dataset, tiered.",
         "C: 06's list age-coordination test in Radial glia, Neuronal IPC and Glioblast, "
         "restricted to cycle-independent genes, null drawn from cycle-independent genes of "
         "matched expression level."],
        f,
        ["Cluster pseudobulks pool cells of all phases, so gene-phase associations are between "
         "clusters with different phase mixes -- an ecological measure, not per-cell phase "
         "expression.",
         "Phase calls are derived from cell-cycle marker genes; canonical markers' phase lean is "
         "partly circular and serves only as a check.",
         "human_dev has no per-cell phase calls, so it has no G1/S/G2M shares or phase lean.",
         "Age points are donors (5-9 per chemistry); trends are across that many people."],
        ["Export per-cell-phase pseudobulks (cell class x phase) in stage 2 for a direct, "
         "non-ecological phase profile of every gene.",
         "Add an explicit G1 duration estimate (e.g. from S-phase fraction under a steady-state "
         "model) per progenitor type and age."])


def figures(out: C.Output, traj: pd.DataFrame, prof: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None:
        return
    for ds in C.DATASETS:
        rows = []
        for chem in C.CHEMISTRIES:
            n = C.ns(ds, chem)
            tr = C.csv(n, "03_cellcycle/proliferation_trajectory.csv")
            tr = tr[(tr.n_cells >= MIN_POINT_CELLS) & tr.cell_class.isin(PROGENITORS)]
            tr["chemistry"] = chem
            rows.append(tr)
        tr = pd.concat(rows)
        metric = "fraction_cycling" if tr["fraction_cycling"].notna().any() else "cell_cycle_score_mean"
        classes = [c for c in PROGENITORS if c in set(tr.cell_class)]
        if not classes:
            continue
        fig, axes = plt.subplots(1, len(classes), figsize=(3.6 * len(classes), 3), squeeze=False)
        for ax, cls in zip(axes.flat, classes):
            for chem, marker in zip(C.CHEMISTRIES, ["o", "^"]):
                h = tr[(tr.cell_class == cls) & (tr.chemistry == chem)].sort_values("age_pcw")
                ax.plot(h.age_pcw, h[metric], marker=marker, lw=1, label=chem)
            ax.set_title(cls, fontsize=9)
            ax.set_xlabel("age (pcw)", fontsize=8)
            ax.tick_params(labelsize=7)
        axes.flat[0].set_ylabel(metric.replace("_", " "), fontsize=8)
        axes.flat[0].legend(fontsize=7)
        fig.suptitle(f"{ds}: proliferation in progenitors over development", fontsize=9)
        fig.tight_layout()
        out.figure(fig, f"proliferation_{ds}", f"Progenitor proliferation vs age, {ds}")
        plt.close(fig)
    if prof.empty or "rho_s_vs_g2m_v2" not in prof:
        return
    p = prof[(prof.dataset == "cortex")].dropna(subset=["rho_proliferation_v2", "rho_s_vs_g2m_v2"])
    if p.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 5))
    for r in p.itertuples():
        kind = r.gene_set.split(":")[0]
        color = {"list": "#c05621", "seed": "#2b6cb0", "module": "#718096"}.get(kind, "k")
        x = (r.rho_proliferation_v2 + r.rho_proliferation_v3) / 2
        y = (r.rho_s_vs_g2m_v2 + r.rho_s_vs_g2m_v3) / 2
        ax.scatter(x, y, s=18 + 2 * np.sqrt(r.n_genes), color=color, alpha=0.8, edgecolor="none")
        if kind != "module" or abs(y) > 0.15:
            ax.annotate(r.gene_set.split(":")[1], (x, y), fontsize=6, xytext=(3, 2),
                        textcoords="offset points")
    ax.axhline(0, color="grey", lw=0.6)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("mean association with proliferation (Spearman, mean of v2/v3)", fontsize=8)
    ax.set_ylabel("mean S (+) vs G2/M (-) lean (partial Spearman)", fontsize=8)
    ax.set_title("cortex: phase profile of gene sets (orange lists, blue seed panels, grey modules)",
                 fontsize=8)
    fig.tight_layout()
    out.figure(fig, "phase_profile_cortex", "Gene-set proliferation association vs S/G2M lean, cortex")
    plt.close(fig)


if __name__ == "__main__":
    main()
