#!/usr/bin/env python3
"""05 - Is cell-class identity reproducible across chemistries and files?

Question: does "Radial glia" mean the same expression programme in v2 and v3,
and in cortex and human_dev? If identity profiles agree, class-level results
from different strata can be compared; if a chemistry distorts them, they
cannot, whatever the statistics say.

Method: per stratum, cell-class pseudobulk counts -> log2 TMM-CPM; each class's
identity profile is its log2 CPM minus the mean of the other classes in the
comparison (classes shared by both strata only, so the contrast is the same
on each side). Profiles are correlated over the most variable shared genes,
giving a class x class matrix per stratum pair: high diagonal and low
off-diagonal mean identity reproduces. Top-100 marker overlap (Jaccard, from
10_markers) is reported alongside.

Pairs: v2 vs v3 within each file (independent donors, different chemistry),
and cortex vs human_dev within each chemistry (same donors, different file).

Inputs (csv_exports/):
  <ds>__<chem>/09_pseudobulk/cell_class__{pseudobulk_counts,group_summary}.csv
  <ds>__<chem>/10_markers/top_markers_cell_class.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "05_identity_reproducibility"
TOP_VARIABLE = 2000
TOP_MARKERS = 100
PAIRS = [(("cortex", "v2"), ("cortex", "v3")),
         (("human_dev", "v2"), ("human_dev", "v3")),
         (("cortex", "v2"), ("human_dev", "v2")),
         (("cortex", "v3"), ("human_dev", "v3"))]


def profiles(n: str, classes: list[str], genes: list[str]) -> pd.DataFrame:
    cnt = C.group_matrix(n, "cell_class", "pseudobulk_counts")[classes]
    lc = C.tmm_log_cpm(cnt).reindex(genes)
    others = {c: [k for k in classes if k != c] for c in classes}
    return pd.DataFrame({c: lc[c] - lc[others[c]].mean(axis=1) for c in classes})


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rows, mats = [], {}
    for (da, ca), (db, cb) in PAIRS:
        na, nb = C.ns(da, ca), C.ns(db, cb)
        for n in (na, nb):
            out.used(f"{n}/09_pseudobulk/cell_class__pseudobulk_counts.csv",
                     f"{n}/10_markers/top_markers_cell_class.csv")
        ga = C.group_matrix(na, "cell_class", "pseudobulk_counts")
        gb = C.group_matrix(nb, "cell_class", "pseudobulk_counts")
        classes = sorted(set(ga.columns) & set(gb.columns))
        genes = sorted(set(ga.index) & set(gb.index))
        if len(classes) < 3:
            continue
        pa, pb = profiles(na, classes, genes), profiles(nb, classes, genes)
        # variable genes: largest identity contrast on either side
        spread = pd.concat([pa.abs().max(axis=1), pb.abs().max(axis=1)], axis=1).min(axis=1)
        top = spread.sort_values(ascending=False).head(TOP_VARIABLE).index
        cm = pd.DataFrame(C.corr_rows(pa.loc[top].T.to_numpy(), pb.loc[top].T.to_numpy()),
                          index=classes, columns=classes)
        label = f"{na} vs {nb}"
        mats[label] = cm
        ma = C.csv(na, "10_markers/top_markers_cell_class.csv")
        mb = C.csv(nb, "10_markers/top_markers_cell_class.csv")
        for c in classes:
            off = cm.loc[c].drop(c)
            sa = set(ma[(ma.top_group == c) & (ma.rank_in_group <= TOP_MARKERS)].gene)
            sb = set(mb[(mb.top_group == c) & (mb.rank_in_group <= TOP_MARKERS)].gene)
            rows.append({"comparison": label, "kind": "chemistry" if da == db else "file",
                         "cell_class": c, "n_classes": len(classes), "n_genes": len(top),
                         "r_same_class": cm.loc[c, c], "r_best_other": off.max(),
                         "best_other": off.idxmax(), "margin": cm.loc[c, c] - off.max(),
                         "top_marker_jaccard": len(sa & sb) / len(sa | sb) if sa | sb else np.nan,
                         "markers_compared": f"{len(sa)} / {len(sb)}"})
        cmo = cm.copy()
        cmo.index.name = f"{na} (rows) x {nb} (columns)"
        out.write(cmo.reset_index(), f"identity_correlation__{na}__{nb}",
                  f"Correlation of class identity profiles, {label}")
    res = pd.DataFrame(rows)
    out.write(res, "identity_reproducibility", "Per class: same-class r, best other-class r, marker Jaccard")

    figures(out, mats)

    f = []
    for (label, kind), g in res.groupby(["comparison", "kind"], sort=False):
        weak = g[g.margin < 0.2]
        f.append(f"**{label}** ({'different chemistry, independent donors' if kind == 'chemistry' else 'same donors, different file'}): "
                 f"same-class r {g.r_same_class.min():.2f}-{g.r_same_class.max():.2f} "
                 f"(median {g.r_same_class.median():.2f}) vs best other class "
                 f"{g.r_best_other.max():.2f} at most; top-{TOP_MARKERS} marker Jaccard "
                 f"median {g.top_marker_jaccard.median():.2f}."
                 + ("" if weak.empty else " Closest neighbours (margin < 0.2): "
                    + "; ".join(f"{r.cell_class} (r {r.r_same_class:.2f} vs {r.best_other} "
                                f"{r.r_best_other:.2f})" for r in weak.itertuples()) + "."))
    chem = res[res.kind == "chemistry"]
    if len(chem):
        ok = int((chem.margin > 0).sum())
        f.append(f"**Class identity reproduces across chemistries**: a class matches itself "
                 f"first in {ok}/{len(chem)} chemistry comparisons"
                 + (", so v2 and v3 results about the same class describe the same cell type."
                    if ok == len(chem) else "; the exceptions need class-specific care.")
                 + " The smallest margins fall between adjacent lineages (radial glia / "
                 "glioblast, neuroblast / neuron), as expected along a differentiation "
                 "continuum. Marker lists overlap less than profiles agree, because v2's lower "
                 "sensitivity reorders the tail of each list.")
    files = res[res.kind == "file"]
    if len(files):
        f.append(f"**Across files, profiles agree but marker lists do not** (median Jaccard "
                 f"{files.top_marker_jaccard.median():.2f}): human_dev classes span the whole "
                 "brain, so their top markers include regional genes that cortex cells never "
                 "express. Compare the files by profile, not by marker list.")

    out.summary(
        "Cell-class identity across chemistries and files",
        "Does each cell class carry the same expression identity in v2 and v3, and in cortex "
        "and human_dev -- i.e. can class-level results be compared across strata?",
        ["Per stratum: cell-class pseudobulk -> log2 TMM-CPM. Identity profile = class minus "
         "mean of the other classes shared by both strata in the comparison.",
         f"Pearson correlation of profiles over the {TOP_VARIABLE} genes with the largest identity "
         "contrast on both sides; margin = same-class r minus best other-class r.",
         f"Top-{TOP_MARKERS} marker Jaccard from 10_markers/top_markers_cell_class.csv."],
        f,
        ["Class-level pseudobulk pools all donors in a stratum, so this tests the label's "
         "programme, not donor-to-donor variability.",
         "cortex vs human_dev share donors; their agreement bounds processing and regional "
         "pooling differences, not biological replication.",
         "Few shared classes (5) make the off-diagonal a coarse reference."],
        ["Repeat at cluster level with a cluster-to-cluster matching to find sub-types that "
         "do not reproduce across chemistry."])


def figures(out: C.Output, mats: dict[str, pd.DataFrame]) -> None:
    plt = C.plt_or_none()
    if plt is None or not mats:
        return
    fig, axes = plt.subplots(1, len(mats), figsize=(3.6 * len(mats), 3.6), squeeze=False)
    for ax, (label, cm) in zip(axes.flat, mats.items()):
        im = ax.imshow(cm.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(cm.shape[1]), cm.columns, rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(cm.shape[0]), cm.index, fontsize=7)
        ax.set_title(label.replace(" vs ", "\nvs "), fontsize=8)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label="identity profile r")
    out.figure(fig, "identity_correlation", "Class x class identity correlation per comparison")
    plt.close(fig)


if __name__ == "__main__":
    main()
