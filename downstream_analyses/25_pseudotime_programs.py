#!/usr/bin/env python3
"""25 - Gene lists along differentiation (pseudotime), separated from developmental age.

Question: when along the radial glia -> IPC -> neuroblast -> neuron path does
each list switch on or off; and once differentiation is held fixed, does the
list still change with the donor's age? An age trend in 03/06 can be a shift
towards more differentiated cells within a class; at matched pseudotime it
cannot.

Method
  Pseudotime from stage-1 script 20 (principal path through class medians in
  the stored latent space; telencephalon only in human_dev); pseudobulks per
  pseudotime bin and per bin x age from script 09.
  A. Trajectories: log2 TMM-CPM per bin (>= 50 cells), genes >= 5 CPM; each
     gene's Spearman with pseudotime; lists / panels: mean rho vs matched random
     sets; marker check (SOX2, PAX6, VIM fall; NEUROD6, STMN2, SNAP25 rise). Each
     list's switch point = pseudotime where its mean Z profile (smoothed over 3
     bins) first crosses halfway between its minimum and maximum, in the
     direction of its trend. v2 x v3 combined, tiered.
  B. Age at matched pseudotime: within each bin with >= 5 age points of >= 30
     cells, each gene's Spearman with age; averaged over bins (weights = age
     points). Lists: mean of that over members vs matched random sets; v2 x v3
     combined, tiered; compared with 06's class-level age coordination.

Inputs:
  csv_exports/<ds>__<chem>/09_pseudobulk/pseudotime_bin{,_x_age}__{pseudobulk_counts,group_summary}.csv
  csv_exports/<ds>__<chem>/20_pseudotime/pseudotime_by_class.csv
  results/06_gene_list_landscape/list_age_coordination_combined.csv (optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "25_pseudotime_programs"
TITLE = "Gene lists along differentiation pseudotime, and age at matched differentiation"
MIN_BIN_CELLS = 50
MIN_AGE_CELLS = 30
MIN_AGES = 5
MIN_CPM = 5.0
N_RANDOM = 2000
MARKERS_DOWN = ["SOX2", "PAX6", "VIM"]
MARKERS_UP = ["NEUROD6", "STMN2", "SNAP25"]
MAX_LISTED = 8


def bin_center(label: str, n_bins: int) -> float:
    return (int(label.replace("pt", "")) + 0.5) / n_bins


def switch_point(profile: np.ndarray, t: np.ndarray, rising: bool) -> float:
    sm = pd.Series(profile).rolling(3, center=True, min_periods=1).mean().to_numpy()
    mid = (sm.min() + sm.max()) / 2
    hit = np.nonzero(sm >= mid if rising else sm <= mid)[0]
    return float(t[hit[0]]) if hit.size else np.nan


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    gene_rows, set_rows, sw_rows, check_rows, age_rows = [], [], [], [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        if not (C.EXPORTS / n / "09_pseudobulk" / "pseudotime_bin__pseudobulk_counts.csv").exists():
            continue
        out.used(f"{n}/09_pseudobulk/pseudotime_bin__pseudobulk_counts.csv")
        sets = C.analysis_gene_sets(ds, groups=("ndd", "cell_cycle"), modules=False, out=out)
        pc = C.EXPORTS / n / "20_pseudotime" / "pseudotime_by_class.csv"
        if pc.exists():
            out.used(f"{n}/20_pseudotime/pseudotime_by_class.csv")
            q = pd.read_csv(pc)
            check_rows.append({"dataset": ds, "chemistry": chem, "path_ordered": bool(q.path_ordered.all()),
                               "class_medians": ", ".join(f"{r.cell_class} {r.q50:.2f}" for r in q.itertuples())})
        cnt = C.group_matrix(n, "pseudotime_bin", "pseudobulk_counts", min_cells=MIN_BIN_CELLS)
        nb = max(int(c.replace("pt", "")) for c in cnt.columns) + 1 if len(cnt.columns) else 0
        nb = max(nb, 20)
        cols = sorted(cnt.columns)
        t = np.array([bin_center(c, nb) for c in cols])
        if len(cols) < MIN_AGES:
            continue
        lc = C.tmm_log_cpm(cnt[cols])
        lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
        tr = C.gene_age_trends(lc, t).rename(columns={"rho": "rho_pseudotime"}).set_index("gene")
        gene_rows.append(tr.reset_index().assign(dataset=ds, chemistry=chem, n_bins=len(cols)))
        for gname, sign in [(g, -1) for g in MARKERS_DOWN] + [(g, 1) for g in MARKERS_UP]:
            if gname in tr.index:
                check_rows.append({"dataset": ds, "chemistry": chem, "marker": gname, "expected_sign": sign,
                                   "rho_pseudotime": float(tr.loc[gname, "rho_pseudotime"])})
        st = C.set_shift_test(tr.rho_pseudotime, tr.mean_log2cpm, sets, rng, N_RANDOM)
        if not st.empty:
            set_rows.append(st.assign(dataset=ds, chemistry=chem, n_bins=len(cols)))
        Z = pd.DataFrame(C.zscore_rows(lc.to_numpy(float)), index=lc.index)
        for name, members in sets.items():
            g = [x for x in members if x in Z.index]
            if len(g) < 5:
                continue
            prof = Z.loc[g].mean(axis=0).to_numpy()
            rising = float(pd.Series(prof).corr(pd.Series(t), method="spearman")) >= 0
            sw_rows.append({"dataset": ds, "chemistry": chem, "gene_set": name, "direction": "rises" if rising else "falls",
                            "switch_point": switch_point(prof, t, rising),
                            "profile": "|".join(f"{v:.2f}" for v in prof)})
        # B: age at matched pseudotime
        p2 = C.EXPORTS / n / "09_pseudobulk" / "pseudotime_bin_x_age__pseudobulk_counts.csv"
        if not p2.exists():
            continue
        out.used(f"{n}/09_pseudobulk/pseudotime_bin_x_age__pseudobulk_counts.csv")
        cnt2 = C.group_matrix(n, "pseudotime_bin_x_age", "pseudobulk_counts", min_cells=MIN_AGE_CELLS)
        by_bin = {}
        for col in cnt2.columns:
            b, age = C.parse_group(col, 2)
            try:
                a = float(age)
            except ValueError:
                continue
            if not C.excluded(ds, age=a):
                by_bin.setdefault(b, []).append((a, col))
        acc, wts = {}, {}
        used_bins = 0
        for b, items in by_bin.items():
            if len(items) < MIN_AGES:
                continue
            items.sort()
            lb = C.tmm_log_cpm(cnt2[[c for _, c in items]])
            lb = lb.loc[lb.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            r = C.gene_age_trends(lb, np.array([a for a, _ in items])).set_index("gene").rho
            for gname, v in r.dropna().items():
                acc[gname] = acc.get(gname, 0.0) + v * len(items)
                wts[gname] = wts.get(gname, 0.0) + len(items)
            used_bins += 1
        if not acc:
            continue
        mean_rho = pd.Series({g: acc[g] / wts[g] for g in acc})
        st = C.set_shift_test(mean_rho, tr.mean_log2cpm.reindex(mean_rho.index).fillna(lc.mean(axis=1).median()),
                              sets, rng, N_RANDOM)
        if not st.empty:
            age_rows.append(st.assign(dataset=ds, chemistry=chem, n_bins_used=used_bins))
    if not gene_rows:
        out.summary(TITLE, "Skipped: no pseudotime pseudobulks yet.",
                    ["Needs 20_pseudotime/ (stage 1) and 09_pseudobulk/pseudotime_bin*.csv (stage 2)."],
                    ["**Not run**: scripts 20 and 09's pseudotime groupings are new; re-run stages 1 and 2, then "
                     "this analysis."], [], [])
        return
    genes = pd.concat(gene_rows, ignore_index=True)
    sets_a = pd.concat(set_rows, ignore_index=True) if set_rows else pd.DataFrame()
    sw = pd.DataFrame(sw_rows)
    checks = pd.DataFrame(check_rows)
    ages = pd.concat(age_rows, ignore_index=True) if age_rows else pd.DataFrame()
    gcomb = C.combine_chemistries(genes, ["dataset", "gene"], effect="rho_pseudotime", weight="n_bins",
                                  labels=("rises along differentiation", "falls along differentiation"),
                                  carry=("rho_pseudotime",))
    scomb = C.combine_chemistries(sets_a, ["dataset", "gene_set"], weight="n_bins",
                                  labels=("rises along differentiation", "falls along differentiation"),
                                  carry=("mean",)) if len(sets_a) else pd.DataFrame()
    acomb = C.combine_chemistries(ages, ["dataset", "gene_set"], labels=("rises with age", "falls with age"),
                                  carry=("mean", "n_bins_used")) if len(ages) else pd.DataFrame()
    swide = sw.pivot_table(index=["dataset", "gene_set"], columns="chemistry", values="switch_point").reset_index() \
        if len(sw) else pd.DataFrame()
    out.write(checks, "pseudotime_checks", "Path ordering from script 20, and marker genes' correlation with pseudotime")
    out.write(genes, "gene_pseudotime_trends_per_stratum", "Per stratum x gene: Spearman with pseudotime")
    out.write(gcomb, "gene_pseudotime_trends_combined", "Per gene: v2 x v3 combined; tier")
    out.write(sets_a, "set_pseudotime_trends_per_stratum", "Per stratum x set: mean rho with pseudotime vs matched null")
    out.write(scomb, "set_pseudotime_trends_combined", "Per set: v2 x v3 combined; tier")
    out.write(sw, "set_switch_points", "Per stratum x set: switch point (half-way crossing) and mean Z profile over bins")
    out.write(ages, "set_age_at_matched_pseudotime_per_stratum", "Per stratum x set: mean within-bin age rho vs matched null")
    out.write(acomb, "set_age_at_matched_pseudotime_combined", "Per set: v2 x v3 combined; tier")
    figures(out, sw)

    f = []
    if len(checks):
        po = checks.dropna(subset=["path_ordered"]) if "path_ordered" in checks else checks.iloc[0:0]
        mk = checks.dropna(subset=["marker"]) if "marker" in checks else checks.iloc[0:0]
        f.append("**Pseudotime sanity**: path ordered RG < IPC < neuroblast < neuron in "
                 + ", ".join(f"{r.dataset} {r.chemistry} ({'yes' if r.path_ordered else 'NO'}: {r.class_medians})"
                             for r in po.itertuples())
                 + ". Markers (expected sign; rho): " + ", ".join(
                     f"{r.dataset} {r.chemistry} {r.marker} {'+' if r.expected_sign > 0 else '-'} {r.rho_pseudotime:+.2f}"
                     for r in mk.itertuples()) + ".")
    if len(scomb):
        rep = scomb[(scomb.tier != "") & scomb.gene_set.str.startswith(("list:", "seed:"))]
        f.append("**Lists along differentiation** (mean rho with pseudotime v2/v3; switch point v2/v3): " + (
            "; ".join(f"{r.gene_set} {r.direction} ({r.mean_v2:+.2f}/{r.mean_v3:+.2f}; switch "
                      + "/".join(f"{x:.2f}" for x in swide[(swide.dataset == r.dataset) & (swide.gene_set == r.gene_set)]
                                 [[c for c in ("v2", "v3") if c in swide]].to_numpy().ravel())
                      + f"; {r.dataset}; {r.tier})" for r in rep.head(MAX_LISTED * 2).itertuples())
            if len(rep) else "none replicated") + ".")
    if len(acomb):
        rep = acomb[(acomb.tier != "") & acomb.gene_set.str.startswith(("list:", "seed:"))]
        f.append("**Age effect at matched differentiation** (mean within-bin rho with age v2/v3): " + (
            "; ".join(f"{r.dataset} {r.gene_set} {r.direction} ({r.mean_v2:+.2f}/{r.mean_v3:+.2f}; {r.tier})"
                      for r in rep.head(MAX_LISTED * 2).itertuples()) if len(rep) else "none replicated")
                 + ". A list that moves with age here does so beyond the shift towards differentiated cells.")
        p06 = C.RESULTS / "06_gene_list_landscape" / "list_age_coordination_combined.csv"
        if p06.exists():
            out.used("results/06_gene_list_landscape/list_age_coordination_combined.csv")
            c06 = pd.read_csv(p06)
            c06 = c06[(c06.tier.fillna("") != "") & c06.cell_class.isin(["Radial glia", "Neuronal IPC", "Neuroblast", "Neuron"])]
            bits = []
            for (ds, gl), g in c06.groupby(["dataset", "gene_list"]):
                a = acomb[(acomb.dataset == ds) & (acomb.gene_set == f"list:{gl}")]
                if len(a):
                    a = a.iloc[0]
                    bits.append(f"{gl} ({ds}; 06: {', '.join(f'{x.direction} in {x.cell_class}' for x in g.itertuples())}) "
                                f"-> at matched pseudotime {a.direction} {a.mean_v2:+.2f}/{a.mean_v3:+.2f} ({a.tier or 'n.s.'})")
            if bits:
                f.append("**06's class-level list age trends, at matched differentiation**: " + "; ".join(bits[:MAX_LISTED]) + ".")
    out.summary(
        TITLE,
        "Where along radial glia -> IPC -> neuroblast -> neuron does each list switch on or off, and does it "
        "still change with donor age once differentiation is held fixed?",
        ["Pseudotime from stage-1 script 20 (principal path through class medians in the stored latent space; "
         "telencephalon only in human_dev); pseudobulks per bin and per bin x age (script 09).",
         f"A: log2 TMM-CPM per bin (>= {MIN_BIN_CELLS} cells), gene Spearman with pseudotime; lists: mean rho vs "
         f"{N_RANDOM:,} matched random sets; switch point = half-way crossing of the smoothed mean Z profile.",
         f"B: within each bin with >= {MIN_AGES} ages of >= {MIN_AGE_CELLS} cells, gene Spearman with age, averaged "
         "over bins; lists vs matched null; v2 x v3 signed Stouffer, BH per dataset, tiered."],
        f,
        ["The path is fixed by the annotated classes; it orders cells but cannot find branches the labels "
         "do not encode (e.g. interneurons in human_dev's telencephalon follow it too).",
         "Equal-width bins hold unequal numbers of cells; sparse bins are dropped.",
         "Age at matched pseudotime still compares donors, so it carries donor differences as 03 does."],
        ["Branch-aware trajectories (e.g. diffusion pseudotime) from the stored neighbour graph."])


def figures(out: C.Output, sw: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or sw.empty:
        return
    for ds in C.DATASETS:
        h = sw[(sw.dataset == ds) & sw.gene_set.str.startswith(("list:", "seed:"))]
        if h.empty:
            continue
        names = list(dict.fromkeys(h.gene_set))[:12]
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
        for ax, chem in zip(axes, C.CHEMISTRIES):
            for nm in names:
                r = h[(h.gene_set == nm) & (h.chemistry == chem)]
                if r.empty:
                    continue
                prof = np.array([float(x) for x in r.iloc[0].profile.split("|")])
                ax.plot(np.linspace(0, 1, len(prof)), prof, lw=1, label=nm.split(":")[1])
            ax.set_title(f"{ds} {chem}", fontsize=8)
            ax.set_xlabel("pseudotime (bins with enough cells)", fontsize=7)
            ax.tick_params(labelsize=6)
        axes[0].set_ylabel("mean Z of list genes", fontsize=7)
        axes[1].legend(fontsize=5, ncol=2)
        fig.tight_layout()
        out.figure(fig, f"pseudotime_profiles_{ds}", f"Gene-set profiles along pseudotime, {ds}")
        plt.close(fig)


if __name__ == "__main__":
    main()
