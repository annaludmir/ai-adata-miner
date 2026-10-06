#!/usr/bin/env python3
"""16 - A transcriptomic age clock: which genes and lists carry the maturation signal?

Question: can a cell type's developmental age be read from its expression,
using a model trained on one set of donors and tested on another? Which
genes carry that signal, and does a gene list predict age better than random
genes of the same expression level?

03 tests genes one at a time; a clock asks the multivariate question and is
validated on people it never saw (v2 and v3 donors are disjoint).

Method
  Samples: class x age pseudobulks (log2 TMM-CPM, TMM within class, as 03),
  classes with >= 4 age points in both chemistries. Genes >= 5 CPM on
  average in both. Each gene is centred within (chemistry, class) and the
  target is age centred within (chemistry, class), so the clock reads age
  within a cell type, not cell-type identity; genes are scaled by their
  training SD.
  Ridge regression (dual form), penalty chosen by leave-one-age-out
  cross-validation inside the training chemistry. Train v2 -> predict v3 and
  train v3 -> predict v2; accuracy = Spearman of predicted vs true age within
  each class, and pooled over classes (class-centred ages). Null for the full
  clock: ages shuffled within class in the training set (200 times).
  Gene weights: coefficient x training SD; genes in the top 200 by |weight| in
  both directions with the same sign are the consensus clock genes.
  List clocks: the same model on a list's genes only (fixed relative
  penalty), mean pooled accuracy over both directions, against 200 random
  gene sets matched on expression decile; one-sided p (better than random),
  BH per dataset.

Inputs:
  csv_exports/<ds>__<chem>/09_pseudobulk/cell_class_x_age__{pseudobulk_counts,detection_fraction,group_summary}.csv
  csv_exports/<ds>__v2/11_panels/panel_coverage.csv; gene lists
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "16_age_clock"
TITLE = "Transcriptomic age clock within cell types, validated across donor sets"
MIN_AGES = 4
MIN_CPM = 5.0
PENALTIES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]     # x mean diagonal of the kernel
LIST_PENALTY = 1.0
N_SHUFFLE = 200
N_RANDOM_SETS = 200
MIN_LIST_GENES = 10
TOP_WEIGHTS = 200
MAX_LISTED = 8


def stratum_samples(n: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(genes x samples log2 CPM, sample table with class and age)."""
    mats, meta = [], []
    for cls, (lc, _, ages, _) in C.class_age_logcpm(n).items():
        if len(ages) < MIN_AGES:
            continue
        mats.append(lc)
        meta += [{"sample": c, "cell_class": cls, "age": a} for c, a in zip(lc.columns, ages)]
    return pd.concat(mats, axis=1), pd.DataFrame(meta)


def centre(X: np.ndarray, groups: np.ndarray) -> np.ndarray:
    X = X.copy()
    for g in np.unique(groups):
        m = groups == g
        X[..., m] -= X[..., m].mean(axis=-1, keepdims=True)
    return X


def ridge_fit(X: np.ndarray, y: np.ndarray, c: float) -> np.ndarray:
    """Dual ridge: X (n x p), returns beta (p). Penalty = c x mean diag(XX')."""
    K = C.dot(X, X.T)
    lam = c * float(np.trace(K)) / len(y)
    alpha = np.linalg.solve(K + lam * np.eye(len(y)), y)
    return C.dot(X.T, alpha)


def accuracy(pred: np.ndarray, age_c: np.ndarray, cls: np.ndarray) -> tuple[float, dict]:
    per = {}
    for k in np.unique(cls):
        m = cls == k
        if m.sum() >= 3:
            per[k] = float(pd.Series(pred[m]).corr(pd.Series(age_c[m]), method="spearman"))
    pooled = float(pd.Series(pred).corr(pd.Series(age_c), method="spearman"))
    return pooled, per


def choose_penalty(X: np.ndarray, y: np.ndarray, ages: np.ndarray, cls: np.ndarray) -> float:
    """Leave-one-age-out CV (all classes at an age held out together)."""
    best, best_c = -np.inf, PENALTIES[0]
    for c in PENALTIES:
        pred = np.full(len(y), np.nan)
        for a in np.unique(ages):
            test = ages == a
            if test.all():
                continue
            pred[test] = C.dot(X[test], ridge_fit(X[~test], y[~test], c))
        ok = np.isfinite(pred)
        score = accuracy(pred[ok], y[ok], cls[ok])[0]
        if score > best:
            best, best_c = score, c
    return best_c


def prepare(ds: str, out: C.Output):
    data = {}
    for chem in C.CHEMISTRIES:
        n = C.ns(ds, chem)
        out.used(f"{n}/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv")
        data[chem] = stratum_samples(n)
    classes = sorted(set(data["v2"][1].cell_class) & set(data["v3"][1].cell_class))
    genes = None
    for chem in C.CHEMISTRIES:
        lc, meta = data[chem]
        keep = meta.cell_class.isin(classes).to_numpy()
        lc = lc.loc[:, keep]
        ok = set(lc.index[lc.mean(axis=1) >= np.log2(MIN_CPM + 1)])
        genes = ok if genes is None else genes & ok
    genes = sorted(genes)
    prepared = {}
    for chem in C.CHEMISTRIES:
        lc, meta = data[chem]
        keep = meta.cell_class.isin(classes).to_numpy()
        meta = meta[keep].reset_index(drop=True)
        X = centre(lc.loc[genes].to_numpy(float)[:, keep], meta.cell_class.to_numpy()).T   # samples x genes
        y = centre(meta.age.to_numpy(float)[None, :], meta.cell_class.to_numpy())[0]
        prepared[chem] = (X, y, meta)
    return prepared, np.array(genes), classes


def run_direction(train, test, cols=None, c=None, rng=None, n_shuffle=0):
    Xa, ya, ma = train
    Xb, yb, mb = test
    if cols is not None:
        Xa, Xb = Xa[:, cols], Xb[:, cols]
    sd = Xa.std(axis=0)
    sd[sd == 0] = 1.0
    Xa, Xb = Xa / sd, Xb / sd
    cls_a, cls_b = ma.cell_class.to_numpy(), mb.cell_class.to_numpy()
    if c is None:
        c = choose_penalty(Xa, ya, ma.age.to_numpy(float), cls_a)
    beta = ridge_fit(Xa, ya, c)
    pooled, per = accuracy(C.dot(Xb, beta), yb, cls_b)
    null = []
    for _ in range(n_shuffle):
        ys = ya.copy()
        for k in np.unique(cls_a):
            m = np.nonzero(cls_a == k)[0]
            ys[m] = ys[rng.permutation(m)]
        null.append(accuracy(C.dot(Xb, ridge_fit(Xa, ys, c)), yb, cls_b)[0])
    return {"c": c, "beta": beta, "sd": sd, "pooled": pooled, "per_class": per, "null": np.array(null)}


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rng = np.random.default_rng(C.SEED)
    acc_rows, weight_rows, list_rows = [], [], []
    for ds in C.DATASETS:
        prepared, genes, classes = prepare(ds, out)
        C.log(f"  {ds}: {len(genes):,} genes, classes {', '.join(classes)}; samples "
              + ", ".join(f"{c} {len(prepared[c][1])}" for c in C.CHEMISTRIES))
        res = {}
        for tr, te in (("v2", "v3"), ("v3", "v2")):
            r = run_direction(prepared[tr], prepared[te], rng=rng, n_shuffle=N_SHUFFLE)
            res[tr] = r
            p = (np.sum(r["null"] >= r["pooled"]) + 1) / (len(r["null"]) + 1)
            acc_rows.append({"dataset": ds, "train": tr, "test": te, "cell_class": "pooled",
                             "spearman_pred_vs_age": r["pooled"], "shuffle_null_mean": float(r["null"].mean()),
                             "shuffle_p": p, "penalty": r["c"], "n_test_samples": len(prepared[te][1])})
            for k, v in r["per_class"].items():
                acc_rows.append({"dataset": ds, "train": tr, "test": te, "cell_class": k,
                                 "spearman_pred_vs_age": v, "penalty": r["c"],
                                 "n_test_samples": int((prepared[te][2].cell_class == k).sum())})
            w = r["beta"] / r["sd"] * prepared[tr][0].std(axis=0)   # coefficient x training SD
            weight_rows.append(pd.DataFrame({"dataset": ds, "train": tr, "gene": genes, "weight": w}))
        # consensus clock genes
        wv2 = weight_rows[-2].set_index("gene").weight
        wv3 = weight_rows[-1].set_index("gene").weight
        top2 = set(wv2.abs().nlargest(TOP_WEIGHTS).index)
        top3 = set(wv3.abs().nlargest(TOP_WEIGHTS).index)
        cons = [g for g in top2 & top3 if np.sign(wv2[g]) == np.sign(wv3[g])]
        for g in cons:
            weight_rows.append(pd.DataFrame({"dataset": [ds], "train": ["consensus"], "gene": [g],
                                             "weight": [(wv2[g] + wv3[g]) / 2]}))
        # list clocks
        raw_level = C.gene_level(ds).reindex(genes).fillna(0).to_numpy()
        bins = C.level_bins(np.log2(raw_level + 1), 10)
        pos = {g: i for i, g in enumerate(genes)}
        sets = C.analysis_gene_sets(ds, groups=("ndd", "cell_cycle"), modules=False, out=out)
        for name, members in sets.items():
            idx = np.array([pos[g] for g in members if g in pos])
            if idx.size < MIN_LIST_GENES:
                continue
            obs = np.mean([run_direction(prepared[a], prepared[b], cols=idx, c=LIST_PENALTY)["pooled"]
                           for a, b in (("v2", "v3"), ("v3", "v2"))])
            rand = C.matched_draws(bins, idx, N_RANDOM_SETS, rng)
            null = np.array([np.mean([run_direction(prepared[a], prepared[b], cols=np.unique(rs),
                                                    c=LIST_PENALTY)["pooled"]
                                      for a, b in (("v2", "v3"), ("v3", "v2"))]) for rs in rand])
            list_rows.append({"dataset": ds, "gene_set": name, "n_genes": int(idx.size),
                              "accuracy": float(obs), "random_mean": float(null.mean()),
                              "random_sd": float(null.std()),
                              "effect_vs_random_sd": (obs - null.mean()) / null.std() if null.std() > 0 else np.nan,
                              "p_better_than_random": (np.sum(null >= obs) + 1) / (len(null) + 1)})
    acc = pd.DataFrame(acc_rows)
    weights = pd.concat(weight_rows, ignore_index=True)
    lists = pd.DataFrame(list_rows)
    if not lists.empty:
        lists["q"] = np.nan
        for _, ix in lists.groupby("dataset").groups.items():
            lists.loc[ix, "q"] = C.bh(lists.loc[ix, "p_better_than_random"])
        lists = lists.sort_values(["dataset", "p_better_than_random"])
    out.write(acc, "clock_accuracy", "Cross-chemistry age prediction: Spearman of predicted vs true age "
              "(pooled over classes with class-centred ages, and per class); shuffle-null p")
    out.write(weights, "clock_gene_weights", "Gene weights (coefficient x training SD) per training "
              "chemistry, and consensus genes (top by |weight| in both with the same sign)")
    out.write(lists, "list_clocks", "Per gene set: mean cross-chemistry accuracy of a clock built on its "
              "genes only vs random sets matched on expression decile")
    figures(out, acc)

    # ---- findings -----------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        a = acc[(acc.dataset == ds) & (acc.cell_class == "pooled")]
        pc = acc[(acc.dataset == ds) & (acc.cell_class != "pooled")]
        if a.empty:
            continue
        per = pc.groupby("cell_class").spearman_pred_vs_age.apply(lambda s: "/".join(f"{x:+.2f}" for x in s))
        f.append(f"**{ds}: age is predictable across donor sets** -- train v2 -> test v3 rho "
                 f"{a.iloc[0].spearman_pred_vs_age:+.2f} (shuffle p = {C.fmt_p(a.iloc[0].shuffle_p)}), v3 -> v2 "
                 f"{a.iloc[1].spearman_pred_vs_age:+.2f} (p = {C.fmt_p(a.iloc[1].shuffle_p)}); per class "
                 "(v2->v3 / v3->v2): " + ", ".join(f"{k} {v}" for k, v in per.items()) + ".")
        cons = weights[(weights.dataset == ds) & (weights.train == "consensus")].sort_values("weight")
        if len(cons):
            f.append(f"**{ds}: consensus clock genes** ({len(cons)} in the top {TOP_WEIGHTS} of both models, "
                     "same sign) -- rising with age: " + ", ".join(cons[cons.weight > 0].sort_values(
                         "weight", ascending=False).gene.head(12))
                     + "; falling: " + ", ".join(cons[cons.weight < 0].gene.head(12)) + ".")
        l = lists[lists.dataset == ds] if not lists.empty else lists
        if len(l):
            good = l[l.q < 0.05]
            f.append(f"**{ds}: gene sets that predict age better than matched random genes** (mean accuracy "
                     "vs random): " + ("; ".join(f"{r.gene_set} {r.accuracy:+.2f} vs {r.random_mean:+.2f} "
                                                 f"(q = {C.fmt_p(r.q)})" for r in good.head(MAX_LISTED).itertuples())
                                       if len(good) else "none at q < 0.05")
                     + ". Weakest relative to random: " + ", ".join(
                         f"{r.gene_set} {r.accuracy:+.2f} vs {r.random_mean:+.2f}"
                         for r in l.sort_values("effect_vs_random_sd").head(3).itertuples()) + ".")
    out.summary(
        TITLE,
        "Can a cell type's developmental age be predicted from expression by a model trained on "
        "other donors; which genes carry the signal; and do gene lists predict age better than "
        "random genes of the same expression?",
        [f"Class x age pseudobulks (log2 TMM-CPM), classes with >= {MIN_AGES} ages in both chemistries, "
         f"genes >= {MIN_CPM:g} CPM; genes and ages centred within chemistry x class; genes scaled by "
         "training SD.",
         "Ridge regression (dual), penalty by leave-one-age-out CV within the training chemistry; "
         "train v2 -> test v3 and back; Spearman of predicted vs true (class-centred) age. Null: "
         f"ages shuffled within class in training ({N_SHUFFLE} times).",
         f"List clocks: the list's genes only, fixed relative penalty {LIST_PENALTY}; mean accuracy of "
         f"both directions vs {N_RANDOM_SETS} random sets matched on expression decile; one-sided p, "
         "BH per dataset."],
        f,
        ["Few training samples (age points x classes, ~20-40 per chemistry): the clock is a "
         "ranking device, not a calibrated age estimate.",
         "The two chemistries cover different age ranges; within-class centring makes accuracy a "
         "rank agreement inside each test chemistry.",
         "human_dev classes pool regions whose sampling changes with age, so its clock can partly "
         "read region."],
        ["B5 pseudotime would separate developmental age from differentiation state within a class."])


def figures(out: C.Output, acc: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or acc.empty:
        return
    h = acc[acc.cell_class != "pooled"].copy()
    h["label"] = h.dataset + " " + h.cell_class
    m = h.pivot_table(index="label", columns="train", values="spearman_pred_vs_age")
    fig, ax = plt.subplots(figsize=(4.5, 0.3 * len(m) + 1.5))
    for j, (col, mk) in enumerate(zip(m.columns, ["o", "^"])):
        ax.scatter(m[col], range(len(m)), marker=mk, label=f"trained on {col}")
    ax.set_yticks(range(len(m)), m.index, fontsize=7)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("Spearman, predicted vs true age (held-out chemistry)", fontsize=8)
    ax.legend(fontsize=7)
    fig.tight_layout()
    out.figure(fig, "age_clock_accuracy", "Cross-chemistry age-clock accuracy per class")
    plt.close(fig)


if __name__ == "__main__":
    main()
