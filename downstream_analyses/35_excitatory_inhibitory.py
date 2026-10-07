#!/usr/bin/env python3
"""35 - Excitatory vs inhibitory neurons: composition over age, expression differences, age trends.

Question: the Neuron and Neuroblast classes pool glutamatergic (excitatory)
and GABAergic / glycinergic (inhibitory) cells. With each cell called
(stage-2 script 23): how does the inhibitory share change with age, within each
region; which genes and gene lists are higher in inhibitory than in excitatory
cells of the same region and age; and do the lists' age trends (06) hold
within each type?

cortex is the EMX1 (dorsal) lineage, so it holds almost no inhibitory neurons;
its calls serve as a check (few should be inhibitory) and its excitatory
neurons as the reference. The comparisons are human_dev's.

Method
  A. Calls per stratum: shares of excitatory / inhibitory / ambiguous /
     unassigned per class; mean scores per call (script 23's marker check).
  B. Inhibitory share = inhibitory / (excitatory + inhibitory) per age point
     (>= 50 assigned cells), per class, in all regions together and within
     each region; Spearman with age, exact permutation; v2 x v3 signed
     Stouffer, BH, tiered.
  C. Inhibitory vs excitatory, matched: pseudobulks per class:type x region x
     age (human_dev; class:type x age in cortex); log2 TMM-CPM over the class's
     groups; each gene's difference inhibitory - excitatory within every
     region x age point holding >= 50 cells of both, averaged over those
     points (so region and age mix cannot drive it). Gene lists and seed panels:
     the same with set means of gene Z-scores, against random sets matched on
     expression decile. v2 x v3 combined, tiered.
  D. Within-type age trends: per class:type, pseudobulks per age (all regions,
     and telencephalon only in human_dev); gene Spearman with age; list age
     coordination vs matched random genes (as 06); v2 x v3 combined, tiered;
     set beside 06's class-level result.

Inputs: csv_exports/<ds>__<chem>/23_neuron_types/*; gene lists; seed panels;
results/06_gene_list_landscape/list_age_coordination_combined.csv (optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "35_excitatory_inhibitory"
TITLE = "Excitatory vs inhibitory neurons: composition, expression differences and age trends"
SUB = "23_neuron_types"
CLASSES = ["Neuron", "Neuroblast"]
MIN_CELLS = 50
MIN_AGES = 5
MIN_CPM = 5.0
N_RANDOM = 2000
TELENCEPHALON = "Telencephalon"
MAX_LISTED = 8
MARKERS = {"SLC17A6", "SLC17A7", "SLC17A8", "NEUROD2", "NEUROD6", "TBR1", "GAD1", "GAD2", "SLC32A1", "SLC6A5",
           "DLX1", "DLX2", "DLX5", "DLX6"}


def load_counts(n: str, grouping: str) -> pd.DataFrame | None:
    p = C.EXPORTS / n / SUB / f"{grouping}__pseudobulk_counts.csv"
    if not p.exists():
        return None
    cnt = pd.read_csv(p)
    cnt = cnt.set_index(cnt.columns[0])
    gs = pd.read_csv(C.EXPORTS / n / SUB / f"{grouping}__group_summary.csv")
    good = set(gs.loc[gs.n_cells >= MIN_CELLS, "group"].astype(str))
    return cnt[[c for c in cnt.columns if c in good]].apply(pd.to_numeric, errors="coerce")


def parse(label: str) -> dict:
    parts = [x.strip() for x in label.split("|")]
    cls, typ = parts[0].split(":")
    out = {"cell_class": cls, "type": typ, "age": float(parts[-1])}
    out["region"] = parts[1] if len(parts) == 3 else "all"
    return out


def part_b(out: C.Output) -> pd.DataFrame:
    rows = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        for fname, by_region in (("neuron_type_by_age", False), ("neuron_type_by_region_x_age", True)):
            p = C.EXPORTS / n / SUB / f"{fname}.csv"
            if not p.exists():
                continue
            out.used(f"{n}/{SUB}/{fname}.csv")
            t = pd.read_csv(p)
            t["assigned"] = t.excitatory + t.inhibitory
            t = t[(t.assigned >= MIN_CELLS) & ~t.age.map(lambda a: C.excluded(ds, age=a))]
            regcol = [c for c in t.columns if c not in ("age", "cell_class", "n_cells", "excitatory", "inhibitory",
                                                        "ambiguous", "unassigned", "assigned")]
            keys = ["cell_class"] + (regcol[:1] if by_region and regcol else [])
            for k, g in t.groupby(keys):
                k = k if isinstance(k, tuple) else (k,)
                g = g.sort_values("age")
                if g.age.nunique() < MIN_AGES:
                    continue
                share = (g.inhibitory / g.assigned).to_numpy()
                rho = float(C.spearman_rows(share[None, :], g.age.to_numpy())[0])
                pp, _ = C.spearman_perm_p(np.array([rho]), g.age.to_numpy())
                rows.append({"dataset": ds, "chemistry": chem, "cell_class": k[0],
                             "region": k[1] if len(k) > 1 else "all", "n_ages": len(g), "rho_vs_age": rho,
                             "perm_p": float(pp[0]), "share_youngest": float(share[0]), "share_oldest": float(share[-1]),
                             "mean_share": float(g.inhibitory.sum() / g.assigned.sum())})
    return pd.DataFrame(rows)


def part_c(out: C.Output, sets_by_ds: dict, rng: np.random.Generator):
    genes_rows, set_rows = [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        grouping = "neuron_type_x_region_x_age"
        cnt = load_counts(n, grouping)
        if cnt is None:
            grouping = "neuron_type_x_age"
            cnt = load_counts(n, grouping)
        if cnt is None or cnt.empty:
            continue
        out.used(f"{n}/{SUB}/{grouping}__pseudobulk_counts.csv")
        meta = pd.DataFrame([parse(c) for c in cnt.columns], index=cnt.columns)
        meta = meta[~meta.age.map(lambda a: C.excluded(ds, age=a))]
        for cls, mc in meta.groupby("cell_class"):
            pairs = []
            for (reg, age), g in mc.groupby(["region", "age"]):
                t = dict(zip(g.type, g.index))
                if {"excitatory", "inhibitory"} <= set(t):
                    pairs.append((t["inhibitory"], t["excitatory"]))
            if len(pairs) < 3:
                continue
            lc = C.tmm_log_cpm(cnt[list(mc.index)])
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            inh = [a for a, _ in pairs]
            exc = [b for _, b in pairs]
            diff = (lc[inh].to_numpy() - lc[exc].to_numpy()).mean(axis=1)
            genes_rows.append(pd.DataFrame({"dataset": ds, "chemistry": chem, "cell_class": cls, "gene": lc.index,
                                            "log2_inh_minus_exc": diff, "n_pairs": len(pairs),
                                            "mean_log2cpm": lc.mean(axis=1).to_numpy()}))
            Z = pd.DataFrame(C.zscore_rows(lc.to_numpy(float)), index=lc.index, columns=lc.columns)
            zi, ze = Z[inh].to_numpy(), Z[exc].to_numpy()
            dz = pd.Series((zi - ze).mean(axis=1), index=lc.index)
            st = C.set_shift_test(dz, lc.mean(axis=1), sets_by_ds[ds], rng, N_RANDOM)
            if not st.empty:
                set_rows.append(st.assign(dataset=ds, chemistry=chem, cell_class=cls, n_pairs=len(pairs)))
    return (pd.concat(genes_rows, ignore_index=True) if genes_rows else pd.DataFrame(),
            pd.concat(set_rows, ignore_index=True) if set_rows else pd.DataFrame())


def part_d(out: C.Output, sets_by_ds: dict, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        for grouping, scope in (("neuron_type_x_age", "all regions"), ("neuron_type_x_region_x_age", TELENCEPHALON)):
            cnt = load_counts(n, grouping)
            if cnt is None or cnt.empty:
                continue
            out.used(f"{n}/{SUB}/{grouping}__pseudobulk_counts.csv")
            meta = pd.DataFrame([parse(c) for c in cnt.columns], index=cnt.columns)
            meta = meta[~meta.age.map(lambda a: C.excluded(ds, age=a))]
            if scope == TELENCEPHALON:
                meta = meta[meta.region == TELENCEPHALON]
            for (cls, typ), g in meta.groupby(["cell_class", "type"]):
                g = g.sort_values("age")
                if g.age.nunique() < MIN_AGES or g.age.duplicated().any():
                    continue
                lc = C.tmm_log_cpm(cnt[list(g.index)])
                lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
                t = C.gene_age_trends(lc, g.age.to_numpy()).set_index("gene")
                st = C.set_shift_test(t.rho, t.mean_log2cpm, sets_by_ds[ds], rng, N_RANDOM)
                if not st.empty:
                    rows.append(st.assign(dataset=ds, chemistry=chem, cell_class=cls, type=typ, scope=scope,
                                          n_ages=len(g)))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    checks = []
    for ds, chem in C.STRATA:
        p = C.EXPORTS / C.ns(ds, chem) / SUB / "neuron_type_marker_check.csv"
        if p.exists():
            out.used(f"{C.ns(ds, chem)}/{SUB}/neuron_type_marker_check.csv")
            checks.append(pd.read_csv(p).assign(dataset=ds, chemistry=chem))
    if not checks:
        out.summary(TITLE, "Skipped: no neuron-type calls yet.",
                    ["Needs 23_neuron_types/ from stage-2 script 23."],
                    ["**Not run**: script 23 is new; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis."],
                    [], [])
        return
    checks = pd.concat(checks, ignore_index=True)
    sets_by_ds = {ds: C.analysis_gene_sets(ds, groups=("ndd", "marker"), modules=False, out=out) for ds in C.DATASETS}
    share = part_b(out)
    scomb = C.combine_chemistries(share, ["dataset", "cell_class", "region"], effect="rho_vs_age", weight="n_ages",
                                  labels=("inhibitory share rises", "inhibitory share falls"),
                                  carry=("share_youngest", "share_oldest", "mean_share")) if len(share) else pd.DataFrame()
    genes, sets_c = part_c(out, sets_by_ds, rng)
    gcomb = pd.DataFrame()
    if len(genes):
        w = genes.pivot_table(index=["dataset", "cell_class", "gene"], columns="chemistry",
                              values="log2_inh_minus_exc").reset_index()
        if {"v2", "v3"} <= set(w.columns):
            w["call"] = np.where((w.v2 >= 1) & (w.v3 >= 1), "inhibitory-high",
                                 np.where((w.v2 <= -1) & (w.v3 <= -1), "excitatory-high", ""))
        gcomb = w
    ccomb = C.combine_chemistries(sets_c, ["dataset", "cell_class", "gene_set"], weight="n_pairs",
                                  labels=("higher in inhibitory", "higher in excitatory"),
                                  carry=("mean", "n_genes")) if len(sets_c) else pd.DataFrame()
    sets_d = part_d(out, sets_by_ds, rng)
    dcomb = C.combine_chemistries(sets_d, ["dataset", "cell_class", "type", "scope", "gene_set"], weight="n_ages",
                                  labels=("rises with age", "falls with age"), carry=("mean",)) \
        if len(sets_d) else pd.DataFrame()
    out.write(checks, "calls", "Per stratum x class x neuron type: cells, mean scores, marker detection (script 23)")
    out.write(share, "inhibitory_share_per_stratum", "Per stratum x class x region: inhibitory share vs age")
    out.write(scomb, "inhibitory_share_combined", "v2 x v3 combined; tier")
    out.write(genes, "gene_inh_vs_exc_per_stratum", "Per stratum x class x gene: log2 inhibitory - excitatory, "
              "averaged over matched region x age points")
    out.write(gcomb, "gene_inh_vs_exc", "Per class x gene: difference per chemistry; call needs |log2| >= 1 in both")
    out.write(sets_c, "set_inh_vs_exc_per_stratum", "Per stratum x class x set: mean Z difference vs matched null")
    out.write(ccomb, "set_inh_vs_exc_combined", "Per class x set: v2 x v3 combined; tier")
    out.write(sets_d, "set_age_trends_by_type_per_stratum", "Per stratum x class x type x scope x set: mean rho with age")
    out.write(dcomb, "set_age_trends_by_type_combined", "v2 x v3 combined; tier")

    f = []
    ck = checks.groupby(["dataset", "cell_class", "neuron_type"]).n_cells.sum().unstack(fill_value=0)
    f.append("**Calls** (cells excitatory / inhibitory / ambiguous / unassigned, both chemistries): " + "; ".join(
        f"{ds} {cls} " + " / ".join(str(int(r.get(t, 0))) for t in ("excitatory", "inhibitory", "ambiguous", "unassigned"))
        for (ds, cls), r in ck.iterrows()) + ". cortex (EMX1 lineage) should hold almost no inhibitory cells.")
    if len(scomb):
        rep = scomb[scomb.tier != ""]
        allr = scomb[scomb.region == "all"]
        f.append("**Inhibitory share** (all regions; youngest -> oldest, v2 / v3): " + "; ".join(
            f"{r.dataset} {r.cell_class} {r.share_youngest_v2:.0%}->{r.share_oldest_v2:.0%} / "
            f"{r.share_youngest_v3:.0%}->{r.share_oldest_v3:.0%} ({r.tier or 'n.s.'})" for r in allr.itertuples()) + ".")
        reg = rep[rep.region != "all"]
        f.append("**Inhibitory share changing with age within a region**: " + ("; ".join(
            f"{r.dataset} {r.cell_class} in {r.region} {r.direction.replace('inhibitory share ', '')} "
            f"({r.share_youngest_v2:.0%}->{r.share_oldest_v2:.0%} / {r.share_youngest_v3:.0%}->{r.share_oldest_v3:.0%}; {r.tier})"
            for r in reg.itertuples()) if len(reg) else "none replicated") + ".")
    if len(gcomb) and "call" in gcomb:
        for (ds, cls), g in gcomb.groupby(["dataset", "cell_class"]):
            g = g[~g.gene.isin(MARKERS)]
            hi = g[g.call == "inhibitory-high"].assign(m=lambda d: (d.v2 + d.v3) / 2).sort_values("m", ascending=False)
            lo = g[g.call == "excitatory-high"].assign(m=lambda d: (d.v2 + d.v3) / 2).sort_values("m")
            f.append(f"**{ds} {cls}: genes higher in inhibitory / excitatory cells of the same region and age** "
                     f"(|log2| >= 1 in both chemistries; markers used for the calls left out): {len(hi)} / {len(lo)}; "
                     "top inhibitory " + ", ".join(hi.gene.head(10)) + "; top excitatory " + ", ".join(lo.gene.head(10)) + ".")
    if len(ccomb):
        rep = ccomb[(ccomb.tier != "") & ccomb.gene_set.str.startswith(("list:", "seed:"))]
        f.append("**Gene lists by neuron type** (mean Z difference inhibitory - excitatory, v2/v3): " + ("; ".join(
            f"{r.dataset} {r.gene_set} {r.direction} in {r.cell_class} ({r.mean_v2:+.2f}/{r.mean_v3:+.2f}; {r.tier})"
            for r in rep.head(MAX_LISTED * 2).itertuples()) if len(rep) else "none replicated") + ".")
    if len(dcomb):
        rep = dcomb[(dcomb.tier != "") & dcomb.gene_set.str.startswith(("list:", "seed:"))]
        f.append("**List age trends within one neuron type** (mean rho v2/v3): " + ("; ".join(
            f"{r.dataset} {r.gene_set} {r.direction} in {r.type} {r.cell_class.lower()}s ({r.scope}; "
            f"{r.mean_v2:+.2f}/{r.mean_v3:+.2f}; {r.tier})" for r in rep.head(MAX_LISTED * 2).itertuples())
            if len(rep) else "none replicated") + ".")
        p06 = C.RESULTS / "06_gene_list_landscape" / "list_age_coordination_combined.csv"
        if p06.exists():
            out.used("results/06_gene_list_landscape/list_age_coordination_combined.csv")
            c06 = pd.read_csv(p06)
            c06 = c06[(c06.tier.fillna("") != "") & c06.cell_class.isin(CLASSES)]
            bits = []
            for r in c06.itertuples():
                w = dcomb[(dcomb.dataset == r.dataset) & (dcomb.cell_class == r.cell_class)
                          & (dcomb.gene_set == f"list:{r.gene_list}") & (dcomb.scope == "all regions")]
                if len(w):
                    bits.append(f"{r.gene_list} {r.direction} in {r.dataset} {r.cell_class} -> " + ", ".join(
                        f"{x.type} {x.tier or 'n.s.'}" for x in w.itertuples()))
            if bits:
                f.append("**06's class-level list age trends, split by neuron type**: " + "; ".join(bits[:MAX_LISTED * 2]) + ".")
    out.summary(
        TITLE,
        "With each neuron and neuroblast called excitatory or inhibitory: how does the inhibitory share change "
        "with age within each region; which genes and lists differ between the two types at the same region and "
        "age; and do list age trends hold within each type?",
        ["Calls from stage-2 script 23: summed CP10K of excitatory (SLC17A6/7/8, NEUROD2/6, TBR1) and inhibitory "
         "(GAD1/2, SLC32A1, SLC6A5, DLX1/2/5/6) genes, log1p; margin 0.5, minimum 1.0.",
         f"B: inhibitory / (excitatory + inhibitory) per age point (>= {MIN_CELLS} assigned cells), all regions and "
         f"per region; Spearman with age (>= {MIN_AGES} ages), exact permutation; v2 x v3 signed Stouffer, tiered.",
         f"C: log2 TMM-CPM over a class's type x region x age pseudobulks; inhibitory - excitatory within each region "
         f"x age point holding >= {MIN_CELLS} cells of both, averaged; sets vs {N_RANDOM:,} random sets matched on "
         "expression decile; v2 x v3 combined, tiered.",
         "D: per class:type, pseudobulk per age (all regions; telencephalon only in human_dev); gene Spearman with "
         "age; list mean rho vs matched null; v2 x v3 combined, tiered."],
        f,
        ["Marker-based calls: immature neurons with little transmitter expression rely on the lineage genes "
         "(NEUROD / TBR1 vs DLX); monoaminergic and cholinergic neurons fall into 'unassigned'.",
         "cortex holds only excitatory (EMX1-lineage) neurons, so inhibitory results are human_dev's, whose regions "
         "were dissected differently by donor; comparisons within region x age remove that mix from C.",
         "Age points are donors (5-11 per chemistry)."],
        ["Sub-divide inhibitory neurons by origin (LHX6 / MGE, NR2F2 / CGE, MEIS2 / LGE) and excitatory neurons by "
         "layer, with the same pattern."])


if __name__ == "__main__":
    main()
