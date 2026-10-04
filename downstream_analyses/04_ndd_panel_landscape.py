#!/usr/bin/env python3
"""04 - Where and when are neurodevelopmental-disorder gene panels active?

Question: which cell classes preferentially express each NDD panel (ASD high
confidence, ID/DD dominant, epilepsy/DEE, chromatin regulators, synaptic and
channel genes), is that preference consistent across donors and both donor
sets -- and are NDD genes over-represented among genes that change with age?

Part A -- cell-class preference. Per stratum, pseudobulk counts per (cell
  class, age point) are TMM-normalised to log2 CPM and each gene is Z-scored
  across all columns. A panel's score in a column is the mean Z of its genes.
  For class c, at every age point, the difference between c and the mean of
  the other classes is taken; T is the mean of that difference over age points
  (age points are donors, so this is donor-replicated). The null is T for
  random gene sets of the same size matched on expression level, because NDD
  genes are long, neuronal and highly expressed and an unmatched null would
  call that "enrichment". v2 and v3 are combined by signed Stouffer.
Part B -- age trends. Panel genes among the replicated rising / falling genes
  of 03, against all genes tested there (hypergeometric).
Part C -- 17_gsea consolidated: panel x group enrichments significant in both
  chemistries.
Part D -- per NDD gene, its most specific cell class (10_markers) in each
  stratum, and whether v2 and v3 agree.

Inputs (csv_exports/):
  <ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,group_summary}.csv
  <ds>__<chem>/11_panels/panel_coverage.csv
  <ds>__<chem>/10_markers/specificity_cell_class.csv
  <ds>__<chem>/17_gsea/gsea_results.csv
  results/03_age_trends_within_cell_class/age_trends_combined.csv (run 03 first)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import _common as C

SLUG = "04_ndd_panel_landscape"
N_RANDOM = 5000
N_BINS = 10            # expression-level bins for matching random sets
MIN_PANEL_GENES = 5
MIN_AGE_POINTS = 3     # age points at which class c and >= 1 other class coexist
MAX_LISTED = 6         # items per bullet in SUMMARY.md; the CSVs hold the rest


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)

    # ---- A: cell-class preference ------------------------------------------
    pref = []
    ndd_by_ns = {}
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv",
                 f"{n}/11_panels/panel_coverage.csv")
        pan = C.panels(n)
        ndd = pan[(pan.panel_group == "ndd") & pan.exported_in_pseudobulk.astype(bool)]
        ndd_by_ns[n] = ndd
        sets = {pname: list(g["gene"]) for pname, g in ndd.groupby("panel")}
        pref.append(C.set_class_preference(n, ds, chem, sets, rng, n_random=N_RANDOM,
                                           n_bins=N_BINS, min_genes=MIN_PANEL_GENES,
                                           min_age_points=MIN_AGE_POINTS))
        C.log(f"  {n}: scored {ndd.panel.nunique()} NDD panels")
    pref = pd.concat(pref, ignore_index=True)
    out.write(pref, "panel_class_preference_per_stratum",
              "Panel x class: mean Z difference vs other classes over age points, matched-null p")

    comb = C.combine_preference(pref)
    out.write(comb, "panel_class_preference_combined",
              "v2 x v3 combined panel preference per class; tier replicated / supported")

    # ---- B: NDD genes among age-trending genes --------------------------------
    tr_path = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_combined.csv"
    # tier and panels are mostly empty strings; read them as text, not guessed types
    tr = pd.read_csv(C.require(tr_path), low_memory=False,
                     dtype={"tier": str, "panels": str})
    tr[["tier", "panels"]] = tr[["tier", "panels"]].fillna("")
    out.used("results/03_age_trends_within_cell_class/age_trends_combined.csv")
    ndd_all = pd.concat(ndd_by_ns.values())[["panel", "gene"]].drop_duplicates()
    # Highly expressed genes are estimated more precisely and so reach
    # significance more often; NDD genes are highly expressed. The plain
    # hypergeometric therefore over-calls, and the test of record draws random
    # sets matched on expression decile (as in part A).
    hyp = []
    for (ds, c), g in tr.groupby(["dataset", "cell_class"]):
        g = g.reset_index(drop=True)
        level = g["mean_log2cpm"].to_numpy()
        bins = np.digitize(level, np.quantile(level, np.linspace(0, 1, N_BINS + 1)[1:-1]))
        pools = [np.nonzero(bins == b)[0] for b in range(N_BINS)]
        pos = {gene: i for i, gene in enumerate(g.gene)}
        for direction in ("up", "down"):
            is_hit = ((g.tier == "replicated") & (g.direction == direction)).to_numpy()
            nh = int(is_hit.sum())
            for pname, pg in ndd_all.groupby("panel"):
                idx = np.array([pos[x] for x in pg.gene if x in pos])
                if idx.size < MIN_PANEL_GENES:
                    continue
                k = int(is_hit[idx].sum())
                rand = np.stack([rng.choice(pools[bins[i]], size=N_RANDOM) for i in idx], axis=1)
                null = is_hit[rand].sum(axis=1)
                p_matched = (np.sum(null >= k) + 1) / (N_RANDOM + 1)
                M, N = len(g), idx.size
                hyp.append({"dataset": ds, "cell_class": c, "direction": direction,
                            "panel": pname, "panel_genes_tested": N, "n_trending": nh,
                            "n_tested": M, "overlap": k,
                            "expected_matched": float(null.mean()),
                            "expected_unmatched": N * nh / M if M else np.nan,
                            "fold_vs_matched": k / null.mean() if null.mean() > 0 else np.nan,
                            "p_matched": p_matched,
                            "p_hypergeom_unmatched": stats.hypergeom.sf(k - 1, M, nh, N) if k > 0 else 1.0,
                            "genes": "|".join(sorted(g.gene.to_numpy()[idx][is_hit[idx]]))})
    hyp = pd.DataFrame(hyp)
    if not hyp.empty:
        hyp["q"] = np.nan
        for _, idx in hyp.groupby("dataset").groups.items():
            hyp.loc[idx, "q"] = C.bh(hyp.loc[idx, "p_matched"])
        hyp = hyp.sort_values(["dataset", "p_matched"])
    out.write(hyp, "panel_age_trend_enrichment",
              "NDD panel genes among replicated rising/falling genes of 03; p vs "
              "expression-matched random sets (hypergeometric shown for comparison)")

    # ---- C: GSEA consolidated ---------------------------------------------------
    # Optional: parts A, B and D do not depend on it, so a missing 17_gsea output
    # (stage 3 unfinished, or 17 skipped a stratum) skips C instead of the analysis.
    gs_rows, gsea_missing = [], []
    for ds in C.DATASETS:
        res = {}
        for chem in C.CHEMISTRIES:
            n = C.ns(ds, chem)
            p = C.EXPORTS / n / "17_gsea" / "gsea_results.csv"
            if not p.exists():
                gsea_missing.append(n)
                continue
            out.used(f"{n}/17_gsea/gsea_results.csv")
            res[chem] = pd.read_csv(p)
        if len(res) < 2:
            continue
        key = ["grouping", "group", "panel_group", "panel"]
        m = res["v2"].merge(res["v3"], on=key, suffixes=("_v2", "_v3"))
        m["significant_both"] = m.significant_v2.astype(bool) & m.significant_v3.astype(bool)
        m.insert(0, "dataset", ds)
        gs_rows.append(m[key + ["dataset", "nes_v2", "fdr_bh_v2", "nes_v3", "fdr_bh_v3",
                                "significant_v2", "significant_v3", "significant_both"]])
    if gsea_missing:
        C.log(f"  WARNING: no 17_gsea/gsea_results.csv for {', '.join(gsea_missing)} -- "
              "part C covers only datasets with both chemistries present")
    gsea = pd.concat(gs_rows, ignore_index=True) if gs_rows else pd.DataFrame()
    if not gsea.empty:
        out.write(gsea[gsea.significant_v2.astype(bool) | gsea.significant_v3.astype(bool)],
                  "gsea_panel_enrichment_by_chemistry",
                  "17_gsea enrichments significant in either chemistry, with the other alongside")

    # ---- D: per-gene specificity agreement -------------------------------------
    spec_rows = []
    for ds in C.DATASETS:
        sp = {}
        for chem in C.CHEMISTRIES:
            n = C.ns(ds, chem)
            out.used(f"{n}/10_markers/specificity_cell_class.csv")
            s = C.csv(n, "10_markers/specificity_cell_class.csv").set_index("gene")
            sp[chem] = s
        genes = sorted(set(ndd_all.gene) & set(sp["v2"].index) & set(sp["v3"].index))
        panel_of = ndd_all.groupby("gene")["panel"].apply(lambda s: "|".join(sorted(set(s))))
        for gname in genes:
            a, b = sp["v2"].loc[gname], sp["v3"].loc[gname]
            spec_rows.append({"dataset": ds, "gene": gname, "panels": panel_of.get(gname, ""),
                              "top_class_v2": a.top_group, "tau_v2": a.tau_specificity,
                              "top_class_v3": b.top_group, "tau_v3": b.tau_specificity,
                              "same_top_class": a.top_group == b.top_group})
    spec = pd.DataFrame(spec_rows)
    out.write(spec, "ndd_gene_top_class",
              "Per NDD gene: most specific cell class and tau in each chemistry")

    figures(out, pref)

    # ---- findings -------------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        g = comb[comb.dataset == ds]
        rep = g[g.tier == "replicated"]
        sup = g[g.tier == "supported"]
        if rep.empty and sup.empty:
            best = g.iloc[0] if len(g) else None
            f.append(f"**{ds}**: no NDD panel shows a cell-class preference beyond matched "
                     "random genes in both donor sets."
                     + (f" Strongest: {best.panel} in {best.cell_class} (effect v2 "
                        f"{best.effect_v2:+.1f}, v3 {best.effect_v3:+.1f} null SDs)." if best is not None else ""))
            continue
        def consistent(r):
            fv2 = r.frac_higher_v2 if r.direction == "enriched" else 1 - r.frac_higher_v2
            fv3 = r.frac_higher_v3 if r.direction == "enriched" else 1 - r.frac_higher_v3
            return min(fv2, fv3)
        for direction in ("enriched", "depleted"):
            d = rep[rep.direction == direction].reindex(
                rep[rep.direction == direction].stouffer_z.abs().sort_values(ascending=False).index)
            if d.empty:
                continue
            txt = "; ".join(f"{r.panel} in {r.cell_class} ({r.effect_v2:+.1f} / {r.effect_v3:+.1f} "
                            f"null SDs; holds at {consistent(r):.0%} of age points or more)"
                            for r in d.head(MAX_LISTED).itertuples())
            more = f" (+{len(d) - MAX_LISTED} more)" if len(d) > MAX_LISTED else ""
            f.append(f"**{ds}: NDD panels {direction} in a cell class, replicated in both donor "
                     f"sets** (effect v2 / v3 vs expression-matched random genes): {txt}{more}.")
        if not sup.empty:
            f.append(f"**{ds}, supported by the combined test only**: "
                     + "; ".join(f"{r.panel} {r.direction} in {r.cell_class}"
                                 for r in sup.head(MAX_LISTED).itertuples())
                     + (f" (+{len(sup) - MAX_LISTED} more)" if len(sup) > MAX_LISTED else "") + ".")
    for ds in C.DATASETS if not hyp.empty else []:
        h = hyp[hyp.dataset == ds]
        sig = h[h.q < 0.05]
        label = "" if ds == "cortex" else " (region-confounded trends, see 03)"
        if sig.empty:
            top = h.iloc[0]
            f.append(f"**{ds}{label}: NDD panels are not over-represented among age-trending "
                     f"genes** once expression level is matched (no panel x class x direction "
                     f"at q < 0.05; strongest {top.panel} {top.direction} in {top.cell_class}: "
                     f"{top.overlap} vs {top.expected_matched:.1f} expected, "
                     f"p = {C.fmt_p(top.p_matched)}).")
        else:
            f.append(f"**{ds}{label}: NDD panels over-represented among replicated age "
                     "trends**: "
                     + "; ".join(f"{r.panel} {r.direction} in {r.cell_class} "
                                 f"({r.overlap} vs {r.expected_matched:.1f} expected from "
                                 f"expression-matched genes, q = {C.fmt_p(r.q)}; "
                                 f"{r.genes.replace('|', ', ')})"
                                 for r in sig.head(MAX_LISTED).itertuples())
                     + (f" (+{len(sig) - MAX_LISTED} more)" if len(sig) > MAX_LISTED else "") + ".")
    if gsea_missing:
        f.append(f"**GSEA (17) results missing for {', '.join(gsea_missing)}**: part C "
                 + ("skipped entirely" if gsea.empty else "covers only the other dataset")
                 + ". Re-run stage 3 (slurm_03_derived.sh), then this analysis.")
    if not gsea.empty:
        both = gsea[gsea.significant_both & (gsea.panel_group == "ndd")]
        f.append("**GSEA (17) NDD enrichments significant in both chemistries**: "
                 + ("; ".join(f"{r.dataset} {r.panel} in {r.group} ({r.grouping})"
                              for r in both.itertuples()) if len(both) else "none")
                 + f". Across all panels: {int(gsea.significant_both.sum())} of "
                 f"{int((gsea.significant_v2.astype(bool) | gsea.significant_v3.astype(bool)).sum())} "
                 "enrichments significant in either chemistry hold in both.")
    if not spec.empty:
        agree = spec.groupby("dataset").same_top_class.mean()
        f.append("**Per-gene cell-class specificity is chemistry-robust for NDD genes**: "
                 "the most specific class agrees between v2 and v3 for "
                 + ", ".join(f"{v:.0%} in {k}" for k, v in agree.items())
                 + " (ndd_gene_top_class.csv).")

    out.summary(
        "NDD gene panels across cell classes and age",
        "Which cell classes preferentially express each neurodevelopmental-disorder panel, "
        "consistently across donors and chemistries -- and are NDD genes enriched among genes "
        "that change with age?",
        ["Panels: the NDD group of the seed panels (panels/README.md), restricted to genes in "
         "the pseudobulk export.",
         "A: log2 TMM-CPM per (class, age point), genes >= 5 CPM, Z-scored across columns; panel "
         "score = mean Z. T = mean over age points of (class - mean of other classes) -- age "
         f"points are donors. Null: {N_RANDOM:,} random sets, one gene drawn from each panel "
         f"gene's expression decile. Effect in null SDs; v2/v3 by signed Stouffer, BH per dataset.",
         "B: overlap with 03's replicated trends; p from expression-matched random sets drawn "
         "from the genes tested in 03 (plain hypergeometric reported alongside -- it over-calls "
         "highly expressed panels).",
         "C: 17_gsea results joined across chemistries. D: 10_markers top class per gene."],
        f,
        ["Seed panels are short hand-picked lists, not SFARI/DDG2P releases; a missing "
         "enrichment can be a panel-coverage problem (11_panels/panel_coverage.csv).",
         "The matched null controls for expression level, not for gene length or "
         "connectivity; long neuronal genes can still look neuron-enriched for that reason.",
         "Gene-set nulls ignore gene-gene correlation, so p-values are optimistic for "
         "co-regulated panels; the replication requirement is the main safeguard.",
         "human_dev classes pool brain regions (see 03)."],
        ["Swap in full SFARI / DDG2P / Epi25 releases via panels/*.csv and re-run 11, 17 and this.",
         "Add gene-length-matched nulls (gene_id_map has Start/End).",
         "Repeat part A at cluster resolution to localise panels to specific neuron types."])


def figures(out: C.Output, pref: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or pref.empty:
        return
    for ds in C.DATASETS:
        g = pref[pref.dataset == ds]
        if g.empty:
            continue
        fig, axes = plt.subplots(1, 2, figsize=(11, 0.45 * g.panel.nunique() + 2.2), squeeze=False)
        for ax, chem in zip(axes.flat, C.CHEMISTRIES):
            h = g[g.chemistry == chem].pivot(index="panel", columns="cell_class",
                                             values="effect_vs_null_sd")
            if h.empty:
                ax.axis("off")
                continue
            im = ax.imshow(h.to_numpy(), cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
            ax.set_xticks(range(h.shape[1]), h.columns, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(h.shape[0]), h.index, fontsize=7)
            ax.set_title(f"{ds} {chem}", fontsize=9)
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8,
                     label="class preference vs matched random genes (null SDs)")
        out.figure(fig, f"panel_preference_{ds}", f"NDD panel x cell class preference, {ds}")
        plt.close(fig)


if __name__ == "__main__":
    main()
