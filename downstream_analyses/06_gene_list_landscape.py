#!/usr/bin/env python3
"""06 - User gene lists: coverage, overlap, cell-class preference, age coordination.

Question: for each gene list in the gene-list folder (config.GENE_LISTS_DIR),
how much of it can these data see, how does it relate to the other lists, in
which cell classes is it preferentially expressed -- and does the list *as a
group* rise or fall with age inside a cell class, consistently in both donor
sets?

The group-level age test matters because single genes rarely survive
correction with 5-11 age points (03), while a list of 50 genes that all drift
the same way is strong evidence even when no member passes alone.

Parts
  A. Coverage: entries matched to symbols (exact / case / Ensembl id), present
     in each file, exported in the pseudobulk, expressed (>= 5 CPM in a class).
  B. Overlap between lists, and with the seed NDD panels (hypergeometric,
     universe = symbols present in human_dev).
  C. Cell-class preference vs expression-matched random genes (as 04 part A).
  D. Age coordination: per stratum x cell class, the mean Spearman rho with age
     of the list's genes (from 03) against random sets matched on expression
     level; v2 and v3 combined by signed Stouffer and tiered.

Inputs:
  gene lists: config.GENE_LISTS_DIR (AIM_GENE_LISTS overrides)
  csv_exports/_cross_dataset/gene_id_map.csv
  csv_exports/<ds>__<chem>/09_pseudobulk/{gene_selection,cell_class__pseudobulk_counts,
                                          cell_class_x_age__pseudobulk_counts}.csv
  csv_exports/<ds>__<chem>/11_panels/panel_coverage.csv
  results/03_age_trends_within_cell_class/age_trends_per_stratum.csv (run 03 first)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import _common as C

SLUG = "06_gene_list_landscape"
TITLE = "Gene lists: coverage, overlap, cell-class preference and age coordination"
N_RANDOM = 5000
N_BINS = 10
MIN_GENES = 5
MAX_LISTED = 6


def coverage(out: C.Output, lists: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for ds in C.DATASETS:
        for name, genes in lists.items():
            m = C.map_genes(genes, ds)
            found = m[m.match != "missing"]
            syms = set(found.symbol)
            row = {"dataset": ds, "gene_list": name, "n_entries": len(genes),
                   "n_matched": len(syms),
                   "n_matched_by_case": int((m.match == "case").sum()),
                   "n_matched_by_ensembl": int((m.match == "ensembl").sum()),
                   "n_missing": int((m.match == "missing").sum()),
                   "missing_examples": "|".join(m.loc[m.match == "missing", "input"].head(10))}
            for chem in C.CHEMISTRIES:
                n = C.ns(ds, chem)
                sel = (C.csv(n, "09_pseudobulk/gene_selection.csv")
                       .drop_duplicates("symbol").set_index("symbol"))
                cnt = C.group_matrix(n, "cell_class", "pseudobulk_counts")
                cpm = cnt.div(cnt.sum(axis=0), axis=1) * 1e6
                exported = syms & set(cnt.index)
                expressed = {g for g in exported if (cpm.loc[g] >= 5).any()}
                below_cut = {g for g in syms - exported
                             if g in sel.index and str(sel.loc[g, "reason"]).startswith("not")}
                row.update({f"n_exported_{chem}": len(exported),
                            f"n_expressed_{chem}": len(expressed),
                            f"n_dropped_below_umi_cut_{chem}": len(below_cut)})
            rows.append(row)
    df = pd.DataFrame(rows)
    out.write(df, "list_coverage",
              "Per list x dataset: entries matched, exported in pseudobulk, expressed")
    return df


def overlaps(out: C.Output, lists: dict[str, list[str]]) -> pd.DataFrame:
    """Pairwise overlap of user lists with each other and with seed NDD panels."""
    # Universe: genes detected anywhere in human_dev (GeneTotalUMIs > 0), plus
    # every set member -- all ~59k annotated symbols would inflate the folds.
    gm = C.csv("_cross_dataset", "gene_id_map.csv", low_memory=False).query("dataset == 'human_dev'")
    detected = set(gm.loc[pd.to_numeric(gm["GeneTotalUMIs"], errors="coerce") > 0, "symbol"].astype(str))
    known = set(gm["symbol"].astype(str))
    sets = {f"list:{k}": set(v) & known for k, v in C.mapped_lists("human_dev").items()
            if k in lists}
    pan = C.panels(C.ns("human_dev", "v2"))
    for pname, g in pan[pan.panel_group == "ndd"].groupby("panel"):
        sets[f"seed:{pname}"] = set(g.gene) & known
    universe = detected.union(*sets.values())
    names = sorted(sets)
    rows = []
    N = len(universe)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if not (a.startswith("list:") or b.startswith("list:")):
                continue
            A, B = sets[a], sets[b]
            k = len(A & B)
            exp = len(A) * len(B) / N if N else np.nan
            rows.append({"set_a": a, "set_b": b, "n_a": len(A), "n_b": len(B),
                         "overlap": k, "jaccard": k / len(A | B) if A | B else np.nan,
                         "expected": exp, "fold": k / exp if exp else np.nan,
                         "p_hypergeom": stats.hypergeom.sf(k - 1, N, len(A), len(B)) if k else 1.0,
                         "shared_genes": "|".join(sorted(A & B)[:40])})
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q"] = C.bh(df["p_hypergeom"])
        df = df.sort_values("p_hypergeom")
    out.write(df, "list_overlaps",
              "Pairwise overlap of gene lists (and seed NDD panels); universe = genes "
              "detected in human_dev plus all set members")
    return df


def age_coordination(out: C.Output, lists: dict[str, list[str]],
                     rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_per_stratum.csv"
    per = pd.read_csv(C.require(path), low_memory=False, dtype={"panels": str})
    out.used("results/03_age_trends_within_cell_class/age_trends_per_stratum.csv")
    rows = []
    for (ds, chem, cls), g in per.groupby(["dataset", "chemistry", "cell_class"], sort=False):
        g = g.reset_index(drop=True)
        rho = g["spearman_rho"].to_numpy(float)
        level = g["mean_log2cpm"].to_numpy(float)
        bins = np.digitize(level, np.quantile(level, np.linspace(0, 1, N_BINS + 1)[1:-1]))
        pools = [np.nonzero(bins == b)[0] for b in range(N_BINS)]
        pos = {gene: i for i, gene in enumerate(g.gene)}
        mapped = C.mapped_lists(ds)
        for name in lists:
            idx = np.array([pos[x] for x in mapped.get(name, []) if x in pos])
            if idx.size < MIN_GENES:
                continue
            obs = float(rho[idx].mean())
            rand = np.stack([rng.choice(pools[bins[i]], size=N_RANDOM) for i in idx], axis=1)
            null = C.set_mean_rows(rho[:, None], rand)[:, 0]
            sd = null.std()
            p = (np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (N_RANDOM + 1)
            rows.append({"dataset": ds, "chemistry": chem, "cell_class": cls, "gene_list": name,
                         "n_genes": int(idx.size), "n_ages": int(g["n_ages"].iloc[0]),
                         "mean_rho": obs, "null_mean": float(null.mean()), "null_sd": float(sd),
                         "effect_vs_null_sd": (obs - null.mean()) / sd if sd > 0 else np.nan,
                         "frac_genes_rising": float((rho[idx] > 0).mean()),
                         "perm_p": p})
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
                     "n_genes_v2": int(a.n_genes), "n_genes_v3": int(b.n_genes),
                     "mean_rho_v2": a.mean_rho, "effect_v2": a.effect_vs_null_sd, "p_v2": a.perm_p,
                     "mean_rho_v3": b.mean_rho, "effect_v3": b.effect_vs_null_sd, "p_v3": b.perm_p,
                     "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
    comb = pd.DataFrame(comb)
    if not comb.empty:
        comb["combined_q"] = np.nan
        for _, ix in comb.groupby("dataset").groups.items():
            comb.loc[ix, "combined_q"] = C.bh(comb.loc[ix, "combined_p"])
        comb["tier"] = C.replication_tier(comb.effect_v2, comb.p_v2, comb.effect_v3, comb.p_v3,
                                          comb.combined_q)
        comb["direction"] = np.where(comb.stouffer_z > 0, "rises with age", "falls with age")
        comb = comb.sort_values(["dataset", "combined_p"])
    out.write(perstr, "list_age_coordination_per_stratum",
              "Per stratum x class x list: mean rho with age of list genes vs matched random sets")
    out.write(comb, "list_age_coordination_combined",
              "v2 x v3 combined list-level age coordination; tier replicated / supported")
    return perstr, comb


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    lists = C.gene_lists()
    src = C.gene_lists_label()
    if not lists:
        C.log(f"  no gene lists found in {src} -- set AIM_GENE_LISTS or copy them there")
        out.summary(TITLE, "Skipped: no gene lists.", [f"Gene-list folder: `{src}`."],
                    [f"**No gene lists found in `{src}`.** Copy them there or set AIM_GENE_LISTS."],
                    [], [])
        return
    C.log(f"  {len(lists)} gene lists from {src}: "
          + ", ".join(f"{k} ({len(v)})" for k, v in lists.items()))
    out.used("_cross_dataset/gene_id_map.csv")
    rng = np.random.default_rng(C.SEED)

    cov = coverage(out, lists)
    ovl = overlaps(out, lists)

    pref = []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/gene_selection.csv",
                 f"{n}/09_pseudobulk/cell_class__pseudobulk_counts.csv",
                 f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv",
                 f"{n}/11_panels/panel_coverage.csv")
        pref.append(C.set_class_preference(n, ds, chem, C.mapped_lists(ds), rng,
                                           n_random=N_RANDOM, n_bins=N_BINS,
                                           min_genes=MIN_GENES, label="gene_list"))
    pref = pd.concat(pref, ignore_index=True)
    out.write(pref, "list_class_preference_per_stratum",
              "List x class: mean Z difference vs other classes over age points, matched-null p")
    pcomb = C.combine_preference(pref, label="gene_list")
    out.write(pcomb, "list_class_preference_combined",
              "v2 x v3 combined list preference per class; tier replicated / supported")

    _, acomb = age_coordination(out, lists, rng)
    figures(out, pcomb, acomb)

    # ---- findings -----------------------------------------------------------
    f = []
    f.append(f"**{len(lists)} gene lists** from `{src}`: "
             + "; ".join(f"{k} ({len(v)})" for k, v in lists.items()) + ".")
    for ds in C.DATASETS:
        c = cov[cov.dataset == ds]
        lowcov = c[c.n_matched < 0.8 * c.n_entries]
        dropped = c[(c.filter(like="n_dropped_below_umi_cut_").sum(axis=1)) > 0]
        f.append(f"**Coverage in {ds}**: "
                 + "; ".join(f"{r.gene_list} {r.n_matched}/{r.n_entries} matched, "
                             f"{r.n_expressed_v2}/{r.n_expressed_v3} expressed (v2/v3)"
                             for r in c.itertuples()) + "."
                 + ("" if lowcov.empty else " Under 80% matched: "
                    + ", ".join(f"{r.gene_list} (e.g. {r.missing_examples.replace('|', ', ')[:80]})"
                                for r in lowcov.itertuples()) + ".")
                 + ("" if dropped.empty else
                    " Some list genes fell below the pseudobulk UMI cut in these exports ("
                    + ", ".join(dropped.gene_list) + "); re-running stage 2 with the gene-list "
                    "folder in place exports them (lists join the panels)."))
    if not ovl.empty:
        strong = ovl[(ovl.q < 0.05) & (ovl.fold >= 2)]
        f.append("**Overlaps**: "
                 + ("; ".join(f"{r.set_a} & {r.set_b}: {r.overlap} shared "
                              f"(Jaccard {r.jaccard:.2f}, {r.fold:.0f}x expected)"
                              for r in strong.head(MAX_LISTED).itertuples())
                    if len(strong) else "no pair overlaps beyond chance at q < 0.05")
                 + ("" if len(strong) <= MAX_LISTED else f" (+{len(strong) - MAX_LISTED} more)")
                 + ".")
    for ds in C.DATASETS:
        g = pcomb[(pcomb.dataset == ds) & (pcomb.tier == "replicated")] if not pcomb.empty else pcomb
        if g.empty:
            f.append(f"**{ds}: no list shows a replicated cell-class preference.**")
            continue
        parts = []
        for name, h in g.groupby("gene_list"):
            up = h[h.direction == "enriched"].sort_values("stouffer_z", ascending=False)
            dn = h[h.direction == "depleted"].sort_values("stouffer_z")
            bits = []
            if len(up):
                bits.append("enriched in " + ", ".join(
                    f"{r.cell_class} ({r.effect_v2:+.1f}/{r.effect_v3:+.1f})" for r in up.itertuples()))
            if len(dn):
                bits.append("depleted in " + ", ".join(r.cell_class for r in dn.itertuples()))
            parts.append(f"{name}: " + "; ".join(bits))
        f.append(f"**{ds}: replicated cell-class preference** (effect in null SDs, v2/v3): "
                 + " | ".join(parts) + ".")
    for ds in C.DATASETS:
        g = acomb[acomb.dataset == ds] if not acomb.empty else acomb
        rep = g[g.tier == "replicated"] if len(g) else g
        label = "" if ds == "cortex" else " (region-confounded, see 03)"
        if rep.empty:
            top = g.iloc[0] if len(g) else None
            f.append(f"**{ds}{label}: no list moves with age as a group in both donor sets.**"
                     + (f" Strongest: {top.gene_list} in {top.cell_class} (mean rho v2 "
                        f"{top.mean_rho_v2:+.2f}, v3 {top.mean_rho_v3:+.2f}; combined p = "
                        f"{C.fmt_p(top.combined_p)})." if top is not None else ""))
        else:
            f.append(f"**{ds}{label}: lists that move with age as a group, replicated** -- "
                     + "; ".join(f"{r.gene_list} {r.direction} in {r.cell_class} (mean rho "
                                 f"{r.mean_rho_v2:+.2f} / {r.mean_rho_v3:+.2f} vs matched null, "
                                 f"q = {C.fmt_p(r.combined_q)})"
                                 for r in rep.head(MAX_LISTED * 2).itertuples())
                     + ("" if len(rep) <= MAX_LISTED * 2 else f" (+{len(rep) - MAX_LISTED * 2} more)")
                     + ".")

    out.summary(
        TITLE,
        "For each user gene list: what do these data see of it, how does it relate to the "
        "other lists, where is it expressed, and does it move with age as a group?",
        [f"Lists read from `{src}` (CSV with a 'gene' column, or its first column; or one "
         "gene per line). Entries matched to symbols exactly, case-insensitively, or via "
         "Ensembl id.",
         "Overlap: hypergeometric against genes detected in human_dev (plus all set "
         "members), BH across pairs.",
         "Cell-class preference: as 04 part A -- log2 TMM-CPM per (class, age point), gene "
         f"Z-scores, class-minus-others difference averaged over age points; {N_RANDOM:,} random "
         "sets matched on expression decile; v2/v3 signed Stouffer, tiered.",
         "Age coordination: mean Spearman rho with age (from 03, TMM log CPM per class x age) "
         "of the list's expressed genes, against random sets matched on expression decile "
         "within the same stratum x class; two-sided p; v2/v3 signed Stouffer, BH per "
         "dataset, tiered (replicated = both chemistries nominal + combined q < 0.05)."],
        f,
        ["Gene-set nulls ignore gene-gene correlation: co-regulated lists get optimistic "
         "p-values. Replication across independent donor sets is the main safeguard.",
         "Age points are donors (5-11 per chemistry), so a list's age trend is a trend across "
         "that many people; human_dev classes pool regions sampled differently by age.",
         "Lists are only as complete as the exports: genes below the pseudobulk UMI cut are "
         "invisible until stage 2 is re-run with the gene-list folder in place."],
        ["Re-run stages 2-4 with the gene-list folder on the cluster so every list gene is "
         "exported.",
         "Weight genes by evidence (e.g. SFARI score) instead of treating lists as flat sets."])


def figures(out: C.Output, pcomb: pd.DataFrame, acomb: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None:
        return
    for ds in C.DATASETS:
        panels = []
        if not pcomb.empty and (pcomb.dataset == ds).any():
            panels.append(("cell-class preference (combined Z)",
                           pcomb[pcomb.dataset == ds].pivot(index="gene_list", columns="cell_class",
                                                           values="stouffer_z")))
        if not acomb.empty and (acomb.dataset == ds).any():
            panels.append(("age coordination (combined Z; + = rises)",
                           acomb[acomb.dataset == ds].pivot(index="gene_list", columns="cell_class",
                                                           values="stouffer_z")))
        if not panels:
            continue
        nrow = max(h.shape[0] for _, h in panels)
        fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 0.45 * nrow + 2.4),
                                 squeeze=False)
        for ax, (title, h) in zip(axes.flat, panels):
            im = ax.imshow(h.to_numpy(float), cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto")
            ax.set_xticks(range(h.shape[1]), h.columns, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(h.shape[0]), h.index, fontsize=7)
            ax.set_title(f"{ds}: {title}", fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        out.figure(fig, f"gene_lists_{ds}", f"Gene-list preference and age coordination, {ds}")
        plt.close(fig)


if __name__ == "__main__":
    main()
