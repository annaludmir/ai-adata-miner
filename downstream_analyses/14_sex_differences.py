#!/usr/bin/env python3
"""14 - Sex differences in NDD gene lists within cell types (age-adjusted).

Question: autism is diagnosed about four times more often in boys. Within a
cell type, are NDD gene lists expressed differently in male and female
donors at the same developmental age?

Data limits, stated first: donor sex is inferred by 01 (Y-chromosome genes,
XIST). Each age point is about one donor, and males are few: 2-4 all-male age
points per chemistry. A per-chemistry test with replication is therefore not
possible; the analysis pools both chemistries of a dataset and adjusts for
chemistry, then asks whether the difference points the same way in each
chemistry on its own. cortex donors are a subset of human_dev's, so the two
datasets do not replicate each other. Treat results as exploratory.

Method
  Per dataset x cell class: log2 TMM-CPM per class x age point (TMM within
  class and chemistry, as 03); age points whose donors all share one inferred
  sex; genes >= 5 CPM. Each gene is regressed on age and chemistry; its
  residuals are scaled by their SD. A set's score per age point = mean scaled
  residual of its genes; statistic = mean over male points - mean over
  female points. Null: random sets matched on expression decile; sex-linked
  genes (Y chromosome, XIST, TSIX) are removed from sets and from the null
  pool. Second null: the sex labels shuffled among age points within each
  chemistry (every relabelling, or 5,000) -- donors differ in ways that are
  not sex (one donor's progenitors may be more mature), and only this null
  asks whether the male/female split is special among splits of the same
  donors. A set must pass both: p = the larger of the two; BH per dataset;
  tier "consistent" = q < 0.05 and the same sign within each chemistry,
  "pooled" = q < 0.05 only.
  Positive control: the same model per gene (Welch t on residuals); Y genes
  and XIST should top the list if the design can see sex at all.

Inputs:
  csv_exports/<ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,detection_fraction,group_summary}.csv
  csv_exports/_cross_dataset/gene_id_map.csv (chromosomes)
  results/01_data_audit/donor_sex_inferred.csv (run 01 first)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import _common as C

SLUG = "14_sex_differences"
TITLE = "Sex differences in gene lists within cell types (age-adjusted, exploratory)"
N_RANDOM = 2000
MIN_GENES = 5
MIN_PER_SEX = 3        # age points per sex, pooled over chemistries
MIN_CPM = 5.0
MAX_LISTED = 8
MAX_LABEL_PERMS = 5000


def label_permutations(sx: np.ndarray, chem: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Boolean (n_perm x points) 'male' masks: sex labels shuffled within each chemistry.

    All within-chemistry relabellings when there are few, else a random sample.
    """
    per_chem = []
    for c in C.CHEMISTRIES:
        idx = np.nonzero(chem == c)[0]
        k = int((sx[idx] == "M").sum())
        per_chem.append([np.array(x, dtype=int) for x in itertools.combinations(idx, k)])
    total = int(np.prod([len(x) for x in per_chem]))
    if total <= MAX_LABEL_PERMS:
        combos = itertools.product(*per_chem)
        masks = np.zeros((total, len(sx)), dtype=bool)
        for i, parts in enumerate(combos):
            for part in parts:
                masks[i, part] = True
        return masks
    masks = np.zeros((MAX_LABEL_PERMS, len(sx)), dtype=bool)
    for i in range(MAX_LABEL_PERMS):
        for opts in per_chem:
            masks[i, opts[rng.integers(len(opts))]] = True
    return masks


def sex_by_age(sex: pd.DataFrame, n: str) -> dict[float, str]:
    s = sex[sex.stratum == n]
    out = {}
    for age, g in s.groupby("age_pcw"):
        calls = set(g.inferred_sex.astype(str))
        if len(calls) == 1 and calls <= {"M", "F"}:
            out[round(float(age), 2)] = calls.pop()
    return out


def class_tables(ds: str, sex: pd.DataFrame, out: C.Output) -> dict:
    """{class: (log2 CPM genes x points, ages, chem, sex)} pooled over chemistries."""
    parts: dict[str, list] = {}
    for chem in C.CHEMISTRIES:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv")
        sx = sex_by_age(sex, n)
        for cls, (lc, _, ages, _) in C.class_age_logcpm(n).items():
            keep = [i for i, a in enumerate(ages) if round(float(a), 2) in sx]
            if not keep:
                continue
            sub = lc.iloc[:, keep]
            sub.columns = [f"{chem}|{c}" for c in sub.columns]
            parts.setdefault(cls, []).append(
                (sub, ages[keep], [chem] * len(keep), [sx[round(float(a), 2)] for a in ages[keep]]))
    tables = {}
    for cls, p in parts.items():
        lc = pd.concat([x[0] for x in p], axis=1, join="inner")
        ages = np.concatenate([x[1] for x in p])
        chem = np.concatenate([x[2] for x in p])
        sx = np.concatenate([x[3] for x in p])
        if (sx == "M").sum() >= MIN_PER_SEX and (sx == "F").sum() >= MIN_PER_SEX:
            tables[cls] = (lc, ages, chem, sx)
    return tables


def residuals(lc: pd.DataFrame, ages: np.ndarray, chem: np.ndarray) -> np.ndarray:
    """Scaled residuals of each gene after age (+ chemistry when both present)."""
    cols = [np.ones_like(ages, dtype=float), ages.astype(float)]
    if len(set(chem)) > 1:
        cols.append((chem == "v3").astype(float))
    D = np.stack(cols, axis=1)
    Y = lc.to_numpy(float)
    beta = np.linalg.lstsq(D, Y.T, rcond=None)[0]
    E = Y - C.dot(D, beta).T
    df = max(len(ages) - D.shape[1], 1)
    sd = np.sqrt((E ** 2).sum(axis=1) / df)
    return np.divide(E, sd[:, None], out=np.zeros_like(E), where=sd[:, None] > 0)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    sex_path = C.RESULTS / "01_data_audit" / "donor_sex_inferred.csv"
    sex = pd.read_csv(C.require(sex_path))
    out.used("results/01_data_audit/donor_sex_inferred.csv", "_cross_dataset/gene_id_map.csv")
    chrom = C.chromosome_map()

    set_rows, gene_rows, design = [], [], []
    for ds in C.DATASETS:
        sets = C.analysis_gene_sets(ds, groups=("ndd",), modules=False, out=out)
        for cls, (lc, ages, chem, sx) in class_tables(ds, sex, out).items():
            lc = lc.loc[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)]
            linked = C.sex_linked(lc.index, chrom)
            E = residuals(lc, ages, chem)
            male, female = sx == "M", sx == "F"
            design.append({"dataset": ds, "cell_class": cls, "n_points": len(sx),
                           "n_male": int(male.sum()), "n_female": int(female.sum()),
                           "male_v2": int((male & (chem == "v2")).sum()), "male_v3": int((male & (chem == "v3")).sum()),
                           "female_v2": int((female & (chem == "v2")).sum()),
                           "female_v3": int((female & (chem == "v3")).sum()),
                           "spearman_sex_vs_age": float(pd.Series(male.astype(float)).corr(
                               pd.Series(ages), method="spearman"))})
            masks = label_permutations(sx, chem, rng)
            design[-1]["n_label_perms"] = len(masks)
            mw = masks / masks.sum(axis=1, keepdims=True)
            fw = ~masks / (~masks).sum(axis=1, keepdims=True)
            # gene level (positive control)
            t, p = stats.ttest_ind(E[:, male], E[:, female], axis=1, equal_var=False)
            g = pd.DataFrame({"dataset": ds, "cell_class": cls, "gene": lc.index, "t_male_vs_female": t,
                              "p": p, "sex_linked": linked,
                              "chromosome": [chrom.get(x, "") for x in lc.index]})
            g["q"] = C.bh(g.p)
            g["rank"] = g.p.rank(method="first").astype(int)
            gene_rows.append(g)
            # set level
            keep = ~linked
            Ek = E[keep]
            pos = {x: i for i, x in enumerate(lc.index[keep])}
            bins = C.level_bins(lc[keep].mean(axis=1), 10)
            for name, genes in sets.items():
                idx = np.array([pos[x] for x in genes if x in pos])
                if idx.size < MIN_GENES:
                    continue
                score = Ek[idx].mean(axis=0)
                d = score[male].mean() - score[female].mean()
                null_scores = C.set_mean_rows(Ek, C.matched_draws(bins, idx, N_RANDOM, rng))
                null = null_scores[:, male].mean(axis=1) - null_scores[:, female].mean(axis=1)
                eff, pp, mu, sd = C.null_effect(d, null)
                dl = C.dot(mw, score) - C.dot(fw, score)
                p_lab = float((np.abs(dl) >= abs(d) - 1e-12).mean())
                within = {}
                for c in C.CHEMISTRIES:
                    m, f_ = male & (chem == c), female & (chem == c)
                    within[c] = score[m].mean() - score[f_].mean() if m.any() and f_.any() else np.nan
                set_rows.append({"dataset": ds, "cell_class": cls, "gene_set": name, "n_genes": int(idx.size),
                                 "n_male": int(male.sum()), "n_female": int(female.sum()),
                                 "male_minus_female": d, "effect_vs_null_sd": eff, "perm_p": pp,
                                 "label_perm_p": p_lab, "n_label_perms": len(masks),
                                 "p_both": max(pp, p_lab),
                                 "difference_v2": within["v2"], "difference_v3": within["v3"]})
    design = pd.DataFrame(design)
    genes = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    res = pd.DataFrame(set_rows)
    if not res.empty:
        res["q"] = np.nan
        for _, ix in res.groupby("dataset").groups.items():
            res.loc[ix, "q"] = C.bh(res.loc[ix, "p_both"])
        same = (np.sign(res.difference_v2) == np.sign(res.male_minus_female)) & \
               (np.sign(res.difference_v3) == np.sign(res.male_minus_female))
        res["tier"] = np.where((res.q < 0.05) & same, "consistent", np.where(res.q < 0.05, "pooled", ""))
        res["direction"] = np.where(res.male_minus_female > 0, "higher in males", "higher in females")
        res = res.sort_values(["dataset", "p_both"])
    out.write(design, "design", "Per dataset x class: age points by inferred sex and chemistry; "
              "Spearman of sex (male = 1) with age (entanglement)")
    out.write(res, "set_sex_differences",
              "Per dataset x class x set: male - female mean scaled residual (age, chemistry adjusted) "
              "vs expression-matched random sets (perm_p) and vs shuffled sex labels (label_perm_p); "
              "within-chemistry differences; tier on the larger p")
    if not genes.empty:
        out.write(genes.sort_values(["dataset", "cell_class", "p"]), "gene_sex_differences",
                  "Per gene x class: Welch t of age/chemistry residuals, male vs female (positive "
                  "control: sex-linked genes should rank first)")

    # ---- findings -----------------------------------------------------------
    f = []
    if design.empty:
        f.append(f"**Not testable**: no cell class has >= {MIN_PER_SEX} all-male and all-female age "
                 "points in any dataset.")
    else:
        f.append("**Design** (age points male/female, pooled over chemistries): "
                 + "; ".join(f"{ds}: " + ", ".join(f"{r.cell_class} {r.n_male}/{r.n_female}"
                                                   for r in g.itertuples())
                             for ds, g in design.groupby("dataset")) + ". Sex vs age: "
                 + ", ".join(f"{ds} rho {g.spearman_sex_vs_age.median():+.2f}" for ds, g in design.groupby("dataset"))
                 + " (median over classes; far from 0 means sex and age are entangled). Possible "
                 "sex relabellings per class: " + ", ".join(
                     f"{r.dataset} {r.cell_class} {r.n_label_perms}" for r in design.itertuples())
                 + " -- with fewer than 20, the label-shuffle p cannot reach 0.05.")
    if not genes.empty:
        ctrl = genes[genes.sex_linked]
        top = ctrl.groupby(["dataset", "cell_class"])["rank"].min()
        f.append("**Positive control** (best rank of a Y gene or XIST among all genes, per class): "
                 + "; ".join(f"{ds} {cls} #{int(r)}" for (ds, cls), r in top.items())
                 + ". Rank 1-3 means the design detects real sex differences.")
        auto = genes[~genes.sex_linked & (genes.q < 0.05)]
        f.append("**Autosomal / X genes differing by sex at q < 0.05**: "
                 + ("; ".join(f"{ds} {cls}: " + ", ".join(g.sort_values("p").gene.head(8))
                              for (ds, cls), g in auto.groupby(["dataset", "cell_class"]))
                    if len(auto) else "none") + ".")
    if not res.empty:
        hits = res[res.tier != ""]
        gene_only = res[(res.perm_p < 0.01) & (res.label_perm_p >= 0.05)]
        f.append("**Gene sets differing by sex** (passing both nulls; effect vs random genes in null SDs; "
                 "label-shuffle p; within-chemistry differences v2/v3): "
                 + ("; ".join(f"{r.dataset} {r.gene_set} {r.direction} in {r.cell_class} ({r.effect_vs_null_sd:+.1f}; "
                              f"p_labels {C.fmt_p(r.label_perm_p)}; {r.difference_v2:+.2f}/{r.difference_v3:+.2f}; {r.tier})"
                              for r in hits.head(MAX_LISTED * 2).itertuples())
                    if len(hits) else "none at q < 0.05")
                 + f". {len(gene_only)} set x class tests beat random genes (p < 0.01) but not shuffled "
                 "sex labels (p >= 0.05): there the donors differ, but not by sex.")
    out.summary(
        TITLE,
        "Within a cell type and at the same age, are NDD gene lists expressed differently in male "
        "and female donors?",
        ["Donor sex from 01 (Y genes, XIST); age points kept when all their donors share one sex.",
         f"Per dataset x class, both chemistries pooled (>= {MIN_PER_SEX} age points per sex): genes "
         f">= {MIN_CPM:g} CPM regressed on age and chemistry; residuals scaled by their SD; set score "
         "= mean scaled residual; male - female difference vs "
         f"{N_RANDOM:,} random sets matched on expression decile; sex-linked genes excluded. BH per "
         "dataset; second null shuffles sex labels among age points within chemistry (all "
         f"relabellings or {MAX_LABEL_PERMS:,}); p = the larger of the two; 'consistent' = q < 0.05 "
         "and same sign within each chemistry.",
         "Gene level: Welch t on residuals; sex-linked genes as positive control."],
        f,
        ["Very few male donors (2-4 age points per chemistry): one male donor's peculiarities can "
         "pass for sex. The random-gene null removes donor-wide shifts, not donor-specific programs.",
         "Pooling chemistries means v2 and v3 are not independent replicates here; the within-"
         "chemistry signs are a weak consistency check.",
         "cortex donors are a subset of human_dev donors."],
        ["A larger donor panel, or per-donor pseudobulks from a dataset with balanced sexes."])


if __name__ == "__main__":
    main()
