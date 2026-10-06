#!/usr/bin/env python3
"""22 - Expression measured directly in each cell-cycle phase.

Question: 09 placed genes and lists on an S-vs-G2/M axis from clusters with
different phase mixes -- an ecological estimate. With pseudobulks of the
cells in each phase, inside one cell class: which genes and lists are higher
in S or in G2/M cells, and in cycling vs non-cycling cells? Does that agree
with 09? And do list age trends in progenitors hold *within* a phase, or are
they a change in the phase mix?

Method
  A. Pseudobulk counts per cell class x phase (09, from this version on) ->
     log2 TMM-CPM within each class across its phases; classes with >= 50
     cells in each of G1, S and G2M. Per gene: log2 S / G2M, and log2 cycling
     (mean of G1, S, G2M) / Non-cycling. Gene calls need |log2 ratio| >= 0.5
     with the same sign in both chemistries. Lists / seed panels: mean log2
     ratio vs random sets matched on expression decile; v2 x v3 combined,
     tiered. Agreement with 09's cluster-level S-vs-G2/M lean.
  B. Pseudobulk counts per class x phase x age: within each progenitor class
     and phase with >= 5 age points of >= 30 cells, list age coordination (mean
     rho of member genes with age vs matched random sets, as 06); compared with
     06's class-level result.

Inputs:
  csv_exports/<ds>__<chem>/09_pseudobulk/cell_class_x_phase{,_x_age}__{pseudobulk_counts,group_summary}.csv
  results/09_cell_cycle_programs/gene_phase_map.csv, results/06_gene_list_landscape/
  list_age_coordination_combined.csv (comparisons, optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "22_phase_resolved_expression"
TITLE = "Expression measured in each cell-cycle phase: S vs G2/M, cycling vs not, and age within a phase"
MIN_PHASE_CELLS = 50
MIN_AGE_CELLS = 30
MIN_AGES = 5
MIN_CPM = 5.0
LFC = 0.5
N_RANDOM = 2000
PROGENITORS = ["Radial glia", "Neuronal IPC", "Glioblast"]
PHASES_B = ["G1", "S", "G2M", "Non-cycling"]
MAX_LISTED = 8


def part_a(out: C.Output, sets_by_ds: dict, rng: np.random.Generator):
    genes, sets = [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        if not (C.EXPORTS / n / "09_pseudobulk" / "cell_class_x_phase__pseudobulk_counts.csv").exists():
            continue
        out.used(f"{n}/09_pseudobulk/cell_class_x_phase__pseudobulk_counts.csv")
        cnt = C.group_matrix(n, "cell_class_x_phase", "pseudobulk_counts", min_cells=MIN_PHASE_CELLS)
        by_class: dict[str, dict[str, str]] = {}
        for col in cnt.columns:
            cls, ph = C.parse_group(col, 2)
            by_class.setdefault(cls, {})[ph] = col
        for cls, cols in sorted(by_class.items()):
            if not {"G1", "S", "G2M"} <= set(cols):
                continue
            lc = C.tmm_log_cpm(cnt[list(cols.values())])
            lc.columns = list(cols.keys())
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            d = pd.DataFrame({"gene": lc.index, "mean_log2cpm": lc.mean(axis=1).to_numpy(),
                              "lfc_s_vs_g2m": (lc["S"] - lc["G2M"]).to_numpy()})
            if "Non-cycling" in lc:
                d["lfc_cycling_vs_noncycling"] = (lc[["G1", "S", "G2M"]].mean(axis=1) - lc["Non-cycling"]).to_numpy()
            d.insert(0, "cell_class", cls)
            d.insert(0, "chemistry", chem)
            d.insert(0, "dataset", ds)
            genes.append(d)
            g = d.set_index("gene")
            for stat in ("lfc_s_vs_g2m", "lfc_cycling_vs_noncycling"):
                if stat not in g:
                    continue
                st = C.set_shift_test(g[stat], g.mean_log2cpm, sets_by_ds[ds], rng, N_RANDOM)
                if not st.empty:
                    sets.append(st.assign(dataset=ds, chemistry=chem, cell_class=cls, statistic=stat))
    return (pd.concat(genes, ignore_index=True) if genes else pd.DataFrame(),
            pd.concat(sets, ignore_index=True) if sets else pd.DataFrame())


def gene_calls(genes: pd.DataFrame) -> pd.DataFrame:
    vals = [c for c in ("lfc_s_vs_g2m", "lfc_cycling_vs_noncycling") if c in genes]
    w = genes.pivot_table(index=["dataset", "cell_class", "gene"], columns="chemistry", values=vals)
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()

    def call(prefix, up, down):
        a, b = w.get(f"{prefix}_v2"), w.get(f"{prefix}_v3")
        if a is None or b is None:
            return ""
        return np.where((a >= LFC) & (b >= LFC), up, np.where((a <= -LFC) & (b <= -LFC), down, ""))

    w["phase_call"] = call("lfc_s_vs_g2m", "S", "G2/M")
    w["cycling_call"] = call("lfc_cycling_vs_noncycling", "cycling", "non-cycling")
    return w


def agreement_with_09(out: C.Output, calls: pd.DataFrame) -> pd.DataFrame:
    p = C.RESULTS / "09_cell_cycle_programs" / "gene_phase_map.csv"
    if not p.exists() or calls.empty or "lfc_s_vs_g2m_v2" not in calls:
        return pd.DataFrame()
    out.used("results/09_cell_cycle_programs/gene_phase_map.csv")
    m = pd.read_csv(p)
    rows = []
    for (ds, cls), g in calls.groupby(["dataset", "cell_class"]):
        mm = m[m.dataset == ds].set_index("gene")
        if "rho_s_vs_g2m_v2" not in mm:
            continue
        g = g.set_index("gene")
        shared = g.index.intersection(mm.index)
        direct = g.loc[shared, ["lfc_s_vs_g2m_v2", "lfc_s_vs_g2m_v3"]].mean(axis=1)
        eco = mm.loc[shared, ["rho_s_vs_g2m_v2", "rho_s_vs_g2m_v3"]].mean(axis=1)
        ok = direct.notna() & eco.notna()
        lean = mm.loc[shared, "phase_lean"].fillna("")
        leaners = lean[lean.isin(["S", "G2/M"])].index
        same = ((lean[leaners] == "S") & (direct[leaners] > 0)) | ((lean[leaners] == "G2/M") & (direct[leaners] < 0))
        rows.append({"dataset": ds, "cell_class": cls, "n_genes": int(ok.sum()),
                     "spearman_direct_vs_09": float(direct[ok].corr(eco[ok], method="spearman")) if ok.sum() > 20 else np.nan,
                     "n_09_leaning_genes": len(leaners), "same_direction_here": int(same.sum())})
    return pd.DataFrame(rows)


def part_b(out: C.Output, sets_by_ds: dict, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        if not (C.EXPORTS / n / "09_pseudobulk" / "cell_class_x_phase_x_age__pseudobulk_counts.csv").exists():
            continue
        out.used(f"{n}/09_pseudobulk/cell_class_x_phase_x_age__pseudobulk_counts.csv")
        cnt = C.group_matrix(n, "cell_class_x_phase_x_age", "pseudobulk_counts", min_cells=MIN_AGE_CELLS)
        parts = {}
        for col in cnt.columns:
            cls, ph, age = C.parse_group(col, 3)
            try:
                a = float(age)
            except ValueError:
                continue
            if cls in PROGENITORS and ph in PHASES_B and not C.excluded(ds, age=a):
                parts.setdefault((cls, ph), []).append((a, col))
        for (cls, ph), items in sorted(parts.items()):
            if len(items) < MIN_AGES:
                continue
            items.sort()
            lc = C.tmm_log_cpm(cnt[[c for _, c in items]])
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            t = C.gene_age_trends(lc, np.array([a for a, _ in items])).set_index("gene")
            st = C.set_shift_test(t.rho, t.mean_log2cpm, sets_by_ds[ds], rng, N_RANDOM)
            if not st.empty:
                rows.append(st.assign(dataset=ds, chemistry=chem, cell_class=cls, phase=ph, n_ages=len(items)))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sets_by_ds = {ds: C.analysis_gene_sets(ds, groups=("ndd", "cell_cycle"), modules=False, out=out)
                  for ds in C.DATASETS}
    genes, sets_a = part_a(out, sets_by_ds, rng)
    if genes.empty:
        out.summary(TITLE, "Skipped: no class x phase pseudobulk yet.",
                    ["Needs 09_pseudobulk/cell_class_x_phase__pseudobulk_counts.csv."],
                    ["**Not run**: script 09 writes the phase groupings from this version of the repo on; "
                     "re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis."], [], [])
        return
    calls = gene_calls(genes)
    labels = {"lfc_s_vs_g2m": ("higher in S", "higher in G2M"),
              "lfc_cycling_vs_noncycling": ("higher in cycling", "higher in non-cycling")}
    comb_a = pd.concat([C.combine_chemistries(g, ["dataset", "cell_class", "statistic", "gene_set"],
                                              labels=labels[stat], carry=("mean",))
                        for stat, g in sets_a.groupby("statistic")], ignore_index=True) if not sets_a.empty else sets_a
    agree = agreement_with_09(out, calls)
    sets_b = part_b(out, sets_by_ds, rng)
    comb_b = C.combine_chemistries(sets_b, ["dataset", "cell_class", "phase", "gene_set"], weight="n_ages",
                                   labels=("rises with age", "falls with age"), carry=("mean",)) \
        if not sets_b.empty else pd.DataFrame()
    out.write(genes, "gene_phase_ratios_per_stratum", "Per stratum x class x gene: log2 S/G2M and cycling/non-cycling")
    out.write(calls, "gene_phase_calls", f"Per class x gene: ratios per chemistry; calls need |log2| >= {LFC} in both")
    out.write(sets_a, "set_phase_ratios_per_stratum", "Per stratum x class x set: mean log2 ratio vs matched null")
    out.write(comb_a, "set_phase_ratios_combined", "Per class x set x statistic: v2 x v3 combined; tier")
    out.write(agree, "agreement_with_09", "Per class: direct S/G2M ratio vs 09's cluster-level lean")
    out.write(sets_b, "set_age_trends_within_phase_per_stratum", "Per stratum x class x phase x set: mean rho with age")
    out.write(comb_b, "set_age_trends_within_phase_combined", "Per class x phase x set: v2 x v3 combined; tier")

    f = []
    cc = calls.groupby(["dataset", "cell_class"])
    f.append("**Genes called by phase** (S / G2M / cycling / non-cycling; same sign, |log2| >= "
             f"{LFC} in both chemistries): " + "; ".join(
                 f"{ds} {cls} {int((g.phase_call == 'S').sum())} / {int((g.phase_call == 'G2/M').sum())} / "
                 f"{int((g.cycling_call == 'cycling').sum())} / {int((g.cycling_call == 'non-cycling').sum())}"
                 for (ds, cls), g in cc) + ".")
    if not agree.empty:
        f.append("**Agreement with 09's cluster-level estimate** (Spearman of direct S/G2M ratio with 09's lean; "
                 "09's S- or G2/M-leaning genes with the same direction here): " + "; ".join(
                     f"{r.dataset} {r.cell_class} rho {r.spearman_direct_vs_09:+.2f}, {r.same_direction_here}/"
                     f"{r.n_09_leaning_genes}" for r in agree.itertuples()) + ".")
    if not comb_a.empty:
        rep = comb_a[(comb_a.tier != "") & ~comb_a.gene_set.str.contains("phase")]
        for stat, label in (("lfc_s_vs_g2m", "S vs G2/M"), ("lfc_cycling_vs_noncycling", "cycling vs non-cycling")):
            r = rep[rep.statistic == stat]
            f.append(f"**Lists and panels by {label}, measured directly** (mean log2 ratio v2/v3): "
                     + ("; ".join(f"{x.gene_set} {x.direction} "
                                  f"in {x.dataset} {x.cell_class} ({x.mean_v2:+.2f}/{x.mean_v3:+.2f}; {x.tier})"
                                  for x in r.head(MAX_LISTED * 2).itertuples()) if len(r) else "none replicated") + ".")
        ctrl = comb_a[comb_a.gene_set.isin(["seed:s_phase", "seed:g2m_phase"]) & (comb_a.statistic == "lfc_s_vs_g2m")]
        if len(ctrl):
            f.append("**Positive control** (seed phase panels, S vs G2/M, mean log2 ratio v2/v3): " + "; ".join(
                f"{x.gene_set} in {x.dataset} {x.cell_class} {x.mean_v2:+.2f}/{x.mean_v3:+.2f}" for x in ctrl.itertuples()) + ".")
    if not comb_b.empty:
        rep = comb_b[comb_b.tier != ""]
        f.append("**List age trends within one phase** (mean rho v2/v3): "
                 + ("; ".join(f"{x.gene_set} {x.direction} in {x.dataset} {x.cell_class} {x.phase} "
                              f"({x.mean_v2:+.2f}/{x.mean_v3:+.2f}; {x.tier})" for x in rep.head(MAX_LISTED * 2).itertuples())
                    if len(rep) else "none replicated") + ".")
        p06 = C.RESULTS / "06_gene_list_landscape" / "list_age_coordination_combined.csv"
        if p06.exists():
            out.used("results/06_gene_list_landscape/list_age_coordination_combined.csv")
            c06 = pd.read_csv(p06)
            c06 = c06[c06.tier.fillna("") != ""]
            held = []
            for x in c06.itertuples():
                w = comb_b[(comb_b.dataset == x.dataset) & (comb_b.cell_class == x.cell_class)
                           & (comb_b.gene_set == f"list:{x.gene_list}")]
                if len(w):
                    same = w[np.sign(w.stouffer_z) == np.sign(x.stouffer_z)]
                    held.append(f"{x.gene_list} in {x.dataset} {x.cell_class}: same direction in "
                                f"{len(same)}/{len(w)} phases, tiered in {int((same.tier != '').sum())}")
            if held:
                f.append("**06's replicated list age trends, phase by phase**: " + "; ".join(held[:MAX_LISTED * 2])
                         + ". A trend that holds within phases is not a change in the phase mix.")
    out.summary(
        TITLE,
        "Inside one cell class, which genes and lists are higher in S or G2/M cells and in cycling vs "
        "non-cycling cells, measured directly; does that agree with 09's cluster-level estimate; and do "
        "list age trends hold within a single phase?",
        [f"Pseudobulk per class x phase (09); log2 TMM-CPM within class; classes with >= {MIN_PHASE_CELLS} cells "
         f"in each of G1, S, G2M; genes >= {MIN_CPM:g} CPM; gene calls |log2| >= {LFC} in both chemistries.",
         f"Lists / seed panels: mean log2 ratio vs {N_RANDOM:,} random sets matched on expression decile; v2 x v3 "
         "signed Stouffer, BH per dataset, tiered.",
         f"Within-phase age trends: class x phase x age pseudobulks, >= {MIN_AGES} ages of >= {MIN_AGE_CELLS} cells; "
         "mean rho of list genes vs matched null; combined and tiered."],
        f,
        ["Phase calls come from cell-cycle marker expression, so phase-panel genes are partly circular "
         "(they define the phases); the informative part is every other gene.",
         "Pseudobulks pool donors; per-chemistry replication is the safeguard."],
        ["Per-donor phase pseudobulks for a donor-level test."])


if __name__ == "__main__":
    main()
