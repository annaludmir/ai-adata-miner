#!/usr/bin/env python3
"""10 - Splicing dynamics: which genes are being switched on or off, where and when?

Question: expression levels say where a gene's mRNA is; the unspliced
(nascent) share says what is happening to it now. A gene being switched on
carries more unspliced RNA than its mature mRNA would predict, one being
switched off carries less (the logic of RNA velocity). So: in which cell
classes along the lineage is each gene, list or module being induced or
shut down -- and inside a cell class, is a gene's unspliced share rising or
falling with age, alongside its expression trend (03)?

The headline use: lists that rise with age in progenitors (06, 09) -- is that
active new transcription (unspliced share rising too), or mRNA accumulating
or stabilising (unspliced share flat or falling)?

Method (cortex only: human_dev has no spliced/unspliced layers)
  * Per group: log unspliced/spliced ratio r = log((U+0.5)/(S+0.5)) from
    pseudobulk counts (script 12), for genes with U+S >= MIN_COUNTS there.
    Comparing a gene with itself across groups cancels its intron content.
  * A gene is compared only across groups where it is clearly on (>= 10 CPM
    and >= 25% of its peak): where a gene is nearly off, its few reads are
    disproportionately unspliced, which made S-phase genes look "switched on"
    in neurons. The remaining coupling (log U/S falls as expression rises, a
    pooled within-gene slope) is subtracted.
  * Each group's median r over its genes is subtracted: a group's overall
    unspliced level shifts for technical reasons (nuclear fraction, stress,
    chemistry), and only relative changes are biology.
  * A. Lineage (cell-class pseudobulk): induction score = a gene's centred r
    in a class minus its mean over classes. Positive = more nascent RNA than
    its mRNA explains (being switched on there); negative = being shut down.
  * B. Age (class x age pseudobulk): Spearman of centred r with age, exact
    permutation p, v2/v3 signed Stouffer, tiered; joined with 03's expression
    trend into: being induced (r up, expression up), induction ahead of
    expression (r up, expression not), being shut down (r down, expression
    down), mRNA catching up (r down, expression up).
  * Gene sets (user lists, seed NDD panels, robust 08 modules): mean score vs
    random sets matched on expression level x baseline unspliced ratio (5 x 5
    bins), v2/v3 combined, BH per part, tiered.

Inputs (csv_exports/):
  cortex__<chem>/12_splicing/{cell_class,cell_class_x_age}__{spliced,unspliced}_counts.csv
  cortex__<chem>/09_pseudobulk/{cell_class,cell_class_x_age}__group_summary.csv
  results/03_age_trends_within_cell_class/age_trends_combined.csv
  gene lists, seed panels, results/08_coexpression_modules/modules.csv (optional)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "10_splicing_dynamics"
TITLE = "Splicing dynamics: where and when genes are being switched on or off"
DATASET = "cortex"
MIN_COUNTS = 20        # U + S in a group for its ratio to be used
MIN_CPM = 10.0         # a gene must be clearly on in a group to be compared there ...
MIN_FRAC_OF_PEAK = 0.25  # ... and at >= this fraction of its peak among the compared groups
MIN_AGES = 5
N_RANDOM = 2000
N_BINS = 5
MAX_LISTED = 8


def ratio_table(n: str, grouping: str):
    """(centred log U/S genes x groups, total counts, pooled raw log U/S per gene).

    Groups need >= MIN_CELLS cells. The pooled raw ratio (all good groups
    summed) is each gene's baseline -- mostly its intron content -- and is
    what the gene-set nulls match on.
    """
    base = C.EXPORTS / n / "12_splicing"
    sp, un = base / f"{grouping}__spliced_counts.csv", base / f"{grouping}__unspliced_counts.csv"
    if not (sp.exists() and un.exists()):
        return None
    S = pd.read_csv(sp, index_col=0)
    U = pd.read_csv(un, index_col=0).reindex(index=S.index, columns=S.columns)
    gs = C.group_summary(n, grouping)
    good = [c for c in S.columns if c in set(gs.loc[gs.n_cells >= C.MIN_CELLS, "group"].astype(str))]
    S, U = S[good].astype(float), U[good].astype(float)
    tot = S + U
    r = np.log((U + 0.5) / (S + 0.5))
    r = r.where(tot >= MIN_COUNTS)
    r = r - r.median(axis=0)                        # remove each group's global level
    pooled = np.log((U.sum(axis=1) + 0.5) / (S.sum(axis=1) + 0.5))
    e = np.log2(tot / tot.sum(axis=0) * 1e6 + 1)    # log2 CPM of U + S per group
    return r, tot, pooled, e


def qualify_and_adjust(r: pd.DataFrame, e: pd.DataFrame, blocks: dict[str, list[str]] | None = None):
    """Keep only groups where a gene is clearly on, and remove the level coupling.

    Where a gene is barely expressed, its few reads are disproportionately
    unspliced, so its unspliced share looks inflated wherever it is nearly off
    (S-phase genes looked "switched on" in neurons). So a gene is compared only
    across groups where it reaches MIN_CPM and MIN_FRAC_OF_PEAK of its peak
    (peak within each block of groups, e.g. one cell class's ages), and the
    remaining coupling -- a pooled within-gene slope of log U/S on log2 CPM --
    is subtracted. Returns (adjusted ratios with NaN where not qualified, slope).
    """
    blocks = blocks or {"all": list(r.columns)}
    R = pd.DataFrame(np.nan, index=r.index, columns=r.columns)
    E = R.copy()
    for cols in blocks.values():
        eb, rb = e[cols], r[cols]
        peak = eb.max(axis=1)
        ok = (eb >= np.log2(MIN_CPM + 1)) & eb.ge(peak + np.log2(MIN_FRAC_OF_PEAK), axis=0) & rb.notna()
        R[cols] = rb.where(ok)
        E[cols] = eb.where(ok)
    num = den = 0.0
    for cols in blocks.values():
        Rb, Eb = R[cols].to_numpy(float), E[cols].to_numpy(float)
        Rd = Rb - np.nanmean(Rb, axis=1, keepdims=True) if np.isfinite(Rb).any() else Rb
        Ed = Eb - np.nanmean(Eb, axis=1, keepdims=True) if np.isfinite(Eb).any() else Eb
        m = np.isfinite(Rd) & np.isfinite(Ed)
        num += float((Rd[m] * Ed[m]).sum())
        den += float((Ed[m] ** 2).sum())
    slope = num / den if den > 0 else 0.0
    return R - slope * E, slope


def matched_null(values: np.ndarray, level: np.ndarray, baseline: np.ndarray,
                 idx: np.ndarray, rng: np.random.Generator) -> np.ndarray | None:
    """Means of random sets matching idx on expression level x baseline ratio bins."""
    ok = np.isfinite(values) & np.isfinite(level) & np.isfinite(baseline)
    if ok.sum() < 50:
        return None
    qb = lambda x: np.digitize(x, np.nanquantile(x[ok], np.linspace(0, 1, N_BINS + 1)[1:-1]))
    bins = np.where(ok, qb(level) * N_BINS + qb(baseline), -1)
    pools = {b: np.nonzero(bins == b)[0] for b in np.unique(bins[idx])}
    if any(p.size == 0 for p in pools.values()):
        return None
    rand = np.stack([rng.choice(pools[bins[i]], size=N_RANDOM) for i in idx], axis=1)
    return C.set_mean_rows(np.nan_to_num(values)[:, None], rand)[:, 0]


def combine(per: pd.DataFrame, keys: list[str], effect: str, weight: str | None) -> pd.DataFrame:
    rows = []
    for k, g in per.groupby(keys, sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        w = [np.sqrt(a[weight]), np.sqrt(b[weight])] if weight else [1.0, 1.0]
        Z, pc = C.signed_stouffer([np.array([a[effect]]), np.array([b[effect]])],
                                  [np.array([a.perm_p]), np.array([b.perm_p])], w)
        rec = dict(zip(keys, k if isinstance(k, tuple) else (k,)))
        rec.update({f"{effect}_v2": a[effect], "p_v2": a.perm_p, f"{effect}_v3": b[effect],
                    "p_v3": b.perm_p, "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
        if "mean" in g:
            rec["mean_v2"], rec["mean_v3"] = a["mean"], b["mean"]
        rows.append(rec)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["combined_q"] = C.bh(df.combined_p)
    df["tier"] = C.replication_tier(df[f"{effect}_v2"], df.p_v2, df[f"{effect}_v3"], df.p_v3,
                                    df.combined_q)
    return df.sort_values("combined_p")


# ---------------------------------------------------------------------------
def lineage(out: C.Output, sets: dict, rng, slopes: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """A. Induction score per gene x class, and per gene set x class."""
    gene_rows, set_rows = [], []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        t = ratio_table(n, "cell_class")
        if t is None:
            continue
        out.used(f"{n}/12_splicing/cell_class__spliced_counts.csv",
                 f"{n}/12_splicing/cell_class__unspliced_counts.csv")
        r, tot, pooled, e = t
        r, slope = qualify_and_adjust(r, e)
        slopes[f"{chem} cell_class"] = slope
        ind = r.sub(r.mean(axis=1), axis=0)
        ind = ind[r.notna().sum(axis=1) >= 2]          # need >= 2 classes where it is clearly on
        level = np.log1p(tot.loc[ind.index].mean(axis=1)).to_numpy()
        baseline = pooled.loc[ind.index].to_numpy()
        long = ind.stack().rename("induction").reset_index()
        long.columns = ["gene", "cell_class", "induction"]
        long.insert(0, "chemistry", chem)
        gene_rows.append(long)
        pos = {g: i for i, g in enumerate(ind.index)}
        for cls in ind.columns:
            v = ind[cls].to_numpy(float)
            for sname, genes in sets.items():
                idx = np.array([pos[x] for x in genes if x in pos and np.isfinite(v[pos[x]])])
                if idx.size < 5:
                    continue
                null = matched_null(v, level, baseline, idx, rng)
                if null is None:
                    continue
                obs = float(v[idx].mean())
                sd = null.std()
                set_rows.append({"chemistry": chem, "cell_class": cls, "gene_set": sname,
                                 "n_genes": int(idx.size), "mean": obs,
                                 "null_mean": float(null.mean()),
                                 "effect": (obs - null.mean()) / sd if sd > 0 else np.nan,
                                 "perm_p": float((np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1)
                                                 / (N_RANDOM + 1))})
    genes = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    if not genes.empty:
        genes = genes.pivot_table(index=["gene", "cell_class"], columns="chemistry",
                                  values="induction").add_prefix("induction_").reset_index()
    out.write(genes, "induction_by_cell_class",
              "Per gene x cell class: centred log U/S minus the gene's mean over classes "
              "(+ = being switched on there), per chemistry")
    per = pd.DataFrame(set_rows)
    comb = combine(per, ["cell_class", "gene_set"], "effect", None) if not per.empty else per
    if not comb.empty:
        comb["direction"] = np.where(comb.stouffer_z > 0, "being switched on", "being switched off")
    out.write(per, "set_induction_by_cell_class_per_stratum",
              "Gene set x class: mean induction score vs matched random sets, per chemistry")
    out.write(comb, "set_induction_by_cell_class",
              "Gene set x class: induction combined over chemistries; tier replicated / supported")
    return genes, comb


def age_dynamics(out: C.Output, sets: dict, rng, slopes: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """B. Unspliced-share trend with age per gene x class, joined with expression trends."""
    gene_rows, set_rows = [], []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        t = ratio_table(n, "cell_class_x_age")
        if t is None:
            continue
        out.used(f"{n}/12_splicing/cell_class_x_age__spliced_counts.csv",
                 f"{n}/12_splicing/cell_class_x_age__unspliced_counts.csv")
        r, tot, pooled, e = t
        by_class = {}
        for col in r.columns:
            cls, age = C.split_class_age(col)
            by_class.setdefault(cls, []).append((age, col))
        # clearly-on and level adjustment within each class's own age series
        r, slope = qualify_and_adjust(r, e, {c: [col for _, col in v] for c, v in by_class.items()})
        slopes[f"{chem} cell_class_x_age"] = slope
        for cls, items in by_class.items():
            items.sort()
            cols = [c for _, c in items]
            ages = np.array([a for a, _ in items])
            if len(ages) < MIN_AGES:
                continue
            sub = r[cols]
            sub = sub[sub.notna().all(axis=1)]           # clearly on at every age point
            if len(sub) < 100:
                continue
            X = sub.to_numpy(float)
            rho = C.spearman_rows(X, ages)
            p, exact = C.spearman_perm_p(rho, ages, n_perm=400_000)
            level = np.log1p(tot.loc[sub.index, cols].mean(axis=1)).to_numpy()
            baseline = pooled.loc[sub.index].to_numpy()
            d = pd.DataFrame({"chemistry": chem, "cell_class": cls, "gene": sub.index,
                              "n_ages": len(ages), "rho": rho, "perm_p": p,
                              "level": level})
            gene_rows.append(d)
            pos = {g: i for i, g in enumerate(sub.index)}
            for sname, genes in sets.items():
                idx = np.array([pos[x] for x in genes if x in pos])
                if idx.size < 5:
                    continue
                null = matched_null(rho, level, baseline, idx, rng)
                if null is None:
                    continue
                obs = float(rho[idx].mean())
                sd = null.std()
                set_rows.append({"chemistry": chem, "cell_class": cls, "gene_set": sname,
                                 "n_genes": int(idx.size), "n_ages": len(ages), "mean": obs,
                                 "null_mean": float(null.mean()),
                                 "effect": (obs - null.mean()) / sd if sd > 0 else np.nan,
                                 "perm_p": float((np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1)
                                                 / (N_RANDOM + 1))})
    per_gene = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    gcomb = pd.DataFrame()
    if not per_gene.empty:
        rows = []
        for (cls, gene), g in per_gene.groupby(["cell_class", "gene"], sort=False):
            if len(g) < 2:
                continue
            g = g.set_index("chemistry")
            a, b = g.loc["v2"], g.loc["v3"]
            rows.append((cls, gene, a.rho, a.perm_p, b.rho, b.perm_p, a.n_ages, b.n_ages))
        gcomb = pd.DataFrame(rows, columns=["cell_class", "gene", "rho_v2", "p_v2", "rho_v3", "p_v3",
                                            "n_ages_v2", "n_ages_v3"])
        Z, pc = C.signed_stouffer([gcomb.rho_v2.to_numpy(), gcomb.rho_v3.to_numpy()],
                                  [gcomb.p_v2.to_numpy(), gcomb.p_v3.to_numpy()],
                                  [np.sqrt(gcomb.n_ages_v2.iloc[0]), np.sqrt(gcomb.n_ages_v3.iloc[0])])
        gcomb["stouffer_z"], gcomb["combined_p"] = Z, pc
        gcomb["combined_q"] = np.nan
        for _, ix in gcomb.groupby("cell_class").groups.items():
            gcomb.loc[ix, "combined_q"] = C.bh(gcomb.loc[ix, "combined_p"])
        gcomb["tier"] = C.replication_tier(gcomb.rho_v2, gcomb.p_v2, gcomb.rho_v3, gcomb.p_v3,
                                           gcomb.combined_q)
        # join expression trends from 03
        tr_path = C.RESULTS / "03_age_trends_within_cell_class" / "age_trends_combined.csv"
        if tr_path.exists():
            out.used("results/03_age_trends_within_cell_class/age_trends_combined.csv")
            tr = pd.read_csv(tr_path, low_memory=False, dtype={"tier": str, "panels": str})
            tr = tr[tr.dataset == DATASET][["cell_class", "gene", "stouffer_z", "tier"]]
            tr = tr.rename(columns={"stouffer_z": "expression_z", "tier": "expression_tier"})
            gcomb = gcomb.merge(tr, on=["cell_class", "gene"], how="left")
            gcomb["expression_tier"] = gcomb["expression_tier"].fillna("")
            up_r = (gcomb.tier == "replicated") & (gcomb.stouffer_z > 0)
            dn_r = (gcomb.tier == "replicated") & (gcomb.stouffer_z < 0)
            ex_up = (gcomb.expression_tier == "replicated") & (gcomb.expression_z > 0)
            ex_dn = (gcomb.expression_tier == "replicated") & (gcomb.expression_z < 0)
            gcomb["dynamics"] = np.select(
                [up_r & ex_up, up_r & ~ex_up & ~ex_dn, dn_r & ex_dn, dn_r & ex_up,
                 up_r & ex_dn, dn_r & ~ex_up & ~ex_dn],
                ["being induced", "induction ahead of expression", "being shut down",
                 "mRNA catching up", "unspliced up, expression down", "shut-down ahead of expression"],
                default="")
        gcomb = gcomb.sort_values(["cell_class", "combined_p"])
    out.write(per_gene, "unspliced_age_trends_per_stratum",
              "Per gene x class x chemistry: Spearman of centred log U/S with age (exact permutation p)")
    out.write(gcomb, "unspliced_age_trends",
              "Per gene x class: unspliced-share age trend combined over chemistries, joined with "
              "03's expression trend into a dynamics category")
    per = pd.DataFrame(set_rows)
    comb = combine(per, ["cell_class", "gene_set"], "effect", "n_ages") if not per.empty else per
    if not comb.empty:
        comb["direction"] = np.where(comb.stouffer_z > 0, "unspliced share rises", "unspliced share falls")
    out.write(per, "set_unspliced_age_trends_per_stratum",
              "Gene set x class: mean unspliced-share age trend vs matched random sets, per chemistry")
    out.write(comb, "set_unspliced_age_trends",
              "Gene set x class: unspliced-share age trend combined over chemistries; tiered")
    return per_gene, gcomb, comb


def global_levels(out: C.Output) -> pd.DataFrame:
    """Each class x age group's overall unspliced level (removed above) -- a QC view."""
    rows = []
    for chem in C.CHEMISTRIES:
        n = C.ns(DATASET, chem)
        base = C.EXPORTS / n / "12_splicing"
        sp, un = base / "cell_class_x_age__spliced_counts.csv", base / "cell_class_x_age__unspliced_counts.csv"
        if not (sp.exists() and un.exists()):
            continue
        S = pd.read_csv(sp, index_col=0)
        U = pd.read_csv(un, index_col=0).reindex(index=S.index, columns=S.columns)
        for col in S.columns:
            s, u = S[col].astype(float), U[col].astype(float)
            ok = (s + u) >= MIN_COUNTS
            if ok.sum() < 100:
                continue
            cls, age = C.split_class_age(col)
            rows.append({"chemistry": chem, "cell_class": cls, "age_pcw": age,
                         "median_log_u_over_s": float(np.median(np.log((u[ok] + 0.5) / (s[ok] + 0.5)))),
                         "overall_unspliced_fraction": float(u.sum() / (u.sum() + s.sum()))})
    df = pd.DataFrame(rows)
    out.write(df, "group_unspliced_levels",
              "Per class x age: overall unspliced level (technical and biological), removed before testing")
    return df


# ---------------------------------------------------------------------------
def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    if not (C.EXPORTS / C.ns(DATASET, "v2") / "12_splicing").exists():
        out.summary(TITLE, "Skipped: no 12_splicing exports for cortex.",
                    ["Needs script 12 (stage 2) outputs."],
                    ["**No splicing exports found** -- run stage 2 (script 12) for cortex."], [], [])
        return
    rng = np.random.default_rng(C.SEED)
    sets = C.analysis_gene_sets(DATASET, groups=("ndd", "cell_cycle"), out=out)

    slopes: dict[str, float] = {}
    genes_a, set_a = lineage(out, sets, rng, slopes)
    _, genes_b, set_b = age_dynamics(out, sets, rng, slopes)
    lev = global_levels(out)
    figures(out, set_a, set_b)

    # ---- findings -----------------------------------------------------------
    f = []
    if not genes_a.empty and {"induction_v2", "induction_v3"} <= set(genes_a.columns):
        agree = genes_a.dropna(subset=["induction_v2", "induction_v3"])
        r = agree.groupby("cell_class").apply(
            lambda d: d.induction_v2.corr(d.induction_v3, method="spearman"), include_groups=False)
        f.append("**Induction scores reproduce across donor sets** (Spearman v2 vs v3 per class): "
                 + ", ".join(f"{k} {v:.2f}" for k, v in r.items())
                 + ". Low values would mean the unspliced signal is noise at this depth.")
    if not set_a.empty:
        rep = set_a[set_a.tier == "replicated"]
        lists = rep[rep.gene_set.str.startswith(("list:", "seed:"))]
        mods = rep[rep.gene_set.str.startswith("module:")]
        f.append("**Where gene sets are being switched on or off along the lineage** (replicated; "
                 "+ = more nascent RNA than their mRNA explains) -- "
                 + ("; ".join(f"{r.gene_set} {r.direction} in {r.cell_class} "
                              f"({r.effect_v2:+.1f}/{r.effect_v3:+.1f} null SDs)"
                              for r in lists.head(MAX_LISTED * 2).itertuples()) if len(lists) else "none")
                 + (f" (+{len(lists) - MAX_LISTED * 2} more)" if len(lists) > MAX_LISTED * 2 else "")
                 + "." + (" Modules: " + "; ".join(f"{r.gene_set.split(':')[1]} {r.direction} in {r.cell_class}"
                                                    for r in mods.head(MAX_LISTED).itertuples()) + "."
                          if len(mods) else ""))
    if not genes_b.empty and "dynamics" in genes_b:
        cnt = genes_b[genes_b.dynamics != ""].groupby(["cell_class", "dynamics"]).size().unstack(fill_value=0)
        n_rep = int((genes_b.tier == "replicated").sum())
        if cnt.empty and n_rep == 0:
            f.append(f"**No single gene's unspliced-share age trend replicates** (0 of "
                     f"{len(genes_b):,} gene x class tests pass in both chemistries at q < 0.05): with "
                     "5-7 age points per chemistry, single genes are underpowered; the gene-set tests "
                     "below pool genes and carry the signal.")
        elif cnt.empty:
            f.append(f"**Genes by dynamics within classes**: {n_rep} gene x class unspliced-share "
                     "trends replicate, but none pairs with a replicated or flat expression trend, so "
                     "none is classed (see `unspliced_age_trends.csv`).")
        else:
            f.append("**Genes by dynamics within classes** (unspliced-share trend and expression "
                     "trend both replicated, or expression flat): "
                     + "; ".join(f"{cls}: " + ", ".join(f"{int(v)} {k}" for k, v in row.items() if v)
                                 for cls, row in cnt.iterrows()) + ".")
    if not set_b.empty:
        rep = set_b[(set_b.tier == "replicated") & set_b.gene_set.str.startswith(("list:", "seed:"))]
        f.append("**Gene sets whose unspliced share changes with age inside a class** (replicated) -- "
                 + ("; ".join(f"{r.gene_set} {r.direction} in {r.cell_class} "
                              f"(mean rho {r.mean_v2:+.2f}/{r.mean_v3:+.2f})"
                              for r in rep.head(MAX_LISTED * 2).itertuples()) if len(rep) else "none")
                 + ". Read with 06/09: a list rising with age whose unspliced share also rises is being "
                 "actively induced; one whose unspliced share is flat is accumulating mRNA.")
    if slopes:
        f.append("**Level coupling removed**: across groups, a gene's log U/S falls as its expression "
                 "rises (pooled within-gene slope per log2 CPM: "
                 + ", ".join(f"{k} {v:+.3f}" for k, v in slopes.items())
                 + "); scores above are net of it, and genes are compared only where clearly on.")
    if not lev.empty:
        trend = lev.groupby(["chemistry", "cell_class"]).apply(
            lambda d: d.median_log_u_over_s.corr(d.age_pcw, method="spearman") if len(d) >= 4 else np.nan,
            include_groups=False).unstack("chemistry")
        f.append("**Overall unspliced level vs age** (removed before every test above; Spearman): "
                 + ", ".join(f"{cls} {row.get('v2', np.nan):+.2f}/{row.get('v3', np.nan):+.2f}"
                             for cls, row in trend.iterrows()) + ".")

    out.summary(
        TITLE,
        "In which cell classes are genes and gene sets being switched on or off (nascent vs mature "
        "RNA), and inside a class, is their unspliced share rising or falling with age alongside "
        "their expression?",
        ["Cortex only (human_dev has no spliced/unspliced layers). Pseudobulk spliced and unspliced "
         f"counts from script 12; groups with >= {C.MIN_CELLS} cells; genes with U+S >= {MIN_COUNTS} "
         "in a group.",
         "Log unspliced/spliced ratio log((U+0.5)/(S+0.5)), each group's median over genes "
         "subtracted (its overall unspliced level shifts for technical reasons); comparing a gene "
         "with itself across groups cancels its intron content.",
         f"A gene is compared only across groups where it is clearly on (>= {MIN_CPM:g} CPM and >= "
         f"{MIN_FRAC_OF_PEAK:.0%} of its peak; for age trends, within the class's own ages): where a "
         "gene is nearly off, its residual reads are disproportionately unspliced. The remaining "
         "coupling between log U/S and expression (a pooled within-gene slope) is subtracted.",
         "A: induction score = centred ratio in a class minus the gene's mean over classes.",
         f"B: Spearman of the centred ratio with age per class (>= {MIN_AGES} ages, complete "
         "profiles), exact permutation p, v2/v3 signed Stouffer, BH per class, tiered; joined with "
         "03's expression trend.",
         f"Gene sets: mean score vs {N_RANDOM:,} random sets matched on expression level x baseline "
         "ratio (5 x 5 bins); v2/v3 combined, BH, tiered."],
        f,
        ["Pseudobulk ratios average over cells in different states; this is a population "
         "measure of induction, not per-cell velocity.",
         "Unspliced counts include intron-retaining and nuclear transcripts; the group-level "
         "centring removes global shifts but not gene-specific changes in retention.",
         "Age points are donors (5-7 per chemistry in cortex).",
         "Cortex only; no replication in human_dev is possible from these files."],
        ["Per-cell velocity (scVelo) on the cortex layers, to place induction along "
         "differentiation pseudotime.",
         "Intron-length-aware matching for the gene-set nulls."])


def figures(out: C.Output, set_a: pd.DataFrame, set_b: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None:
        return
    panels = []
    for df, title in ((set_a, "induction along the lineage"), (set_b, "unspliced-share trend with age")):
        if df is None or df.empty:
            continue
        h = df[df.gene_set.str.startswith(("list:", "seed:"))].pivot(
            index="gene_set", columns="cell_class", values="stouffer_z")
        if not h.empty:
            panels.append((title, h))
    if not panels:
        return
    nrow = max(h.shape[0] for _, h in panels)
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 0.38 * nrow + 2.4), squeeze=False)
    for ax, (title, h) in zip(axes.flat, panels):
        im = ax.imshow(h.to_numpy(float), cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto")
        ax.set_xticks(range(h.shape[1]), h.columns, rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(h.shape[0]), [x.split(":", 1)[1] for x in h.index], fontsize=7)
        ax.set_title(f"cortex: {title} (combined Z)", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    out.figure(fig, "splicing_dynamics_cortex", "Gene-set induction by class and unspliced-share age trends")
    plt.close(fig)


if __name__ == "__main__":
    main()
