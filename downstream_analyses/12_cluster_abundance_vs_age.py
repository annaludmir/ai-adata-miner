#!/usr/bin/env python3
"""12 - Which sub-types (clusters) expand or shrink with age, within their cell class.

Question: 02 found which broad cell classes grow or shrink over development.
Inside a class, do particular sub-types take over -- consistently in both
donor sets? A sub-type's share is taken within its parent class, so a shift
in the class itself (02's result) does not count again here.

Method
  Clusterings whose labels are shared across chemistries (a label's v2
  expression profile best matches the same label in v3): cortex leiden_scVI,
  human_dev cluster_id. (cortex's finer clusterings were computed per chemistry, so a
  label there names different cells in v2 and v3.)
  Cells per cluster x donor = cluster size x the cluster's donor fractions
  (04_clusters). Parent = the cluster's dominant cell class (cortex), or class
  x dominant region (human_dev, whose donors were dissected differently).
  Clusters need >= 50% of cells in that class. Per parent, donors with
  >= 50 cells in it, >= 5 donors: each donor's counts over the parent's
  clusters -> centred log-ratio (composition-aware); Spearman of each
  cluster's CLR with donor age, exact permutation p. v2 x v3: signed
  Stouffer (weights sqrt(donors)), BH per dataset, tiered.

Inputs (csv_exports/):
  <ds>__<chem>/04_clusters/cluster_composition_<clustering>_by_donor.csv
  <ds>__<chem>/04_clusters/cluster_profile_<clustering>.csv
  <ds>__<chem>/05_confounds/crosstab_age_x_donor.csv
  <ds>__<chem>/10_markers/top_markers_cluster_<clustering>.csv (labels only)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

import _common as C

SLUG = "12_cluster_abundance_vs_age"
TITLE = "Sub-type (cluster) abundance over development, within cell classes"
CLUSTERING = {"cortex": "leiden_scVI", "human_dev": "cluster_id"}
MIN_PURITY = 0.5
MIN_PARENT_CELLS = 50
MIN_DONORS = 5
MAX_LISTED = 8


def parent_of(prof: pd.DataFrame, ds: str) -> pd.Series:
    cls = prof["dominant_cell_class"].astype(str)
    return cls + " | " + prof["dominant_region"].astype(str) if ds == "human_dev" else cls


def per_stratum(out: C.Output) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, shares = [], []
    for ds, chem in C.STRATA:
        n, cl = C.ns(ds, chem), CLUSTERING[ds]
        out.used(f"{n}/04_clusters/cluster_composition_{cl}_by_donor.csv",
                 f"{n}/04_clusters/cluster_profile_{cl}.csv", f"{n}/05_confounds/crosstab_age_x_donor.csv")
        comp = C.csv(n, f"04_clusters/cluster_composition_{cl}_by_donor.csv")
        comp["cluster"] = comp["cluster"].astype(str)
        comp = comp.set_index("cluster")
        prof = C.csv(n, f"04_clusters/cluster_profile_{cl}.csv")
        prof["cluster"] = prof["cluster"].astype(str)
        prof = prof.set_index("cluster").reindex(comp.index)
        counts = comp.mul(prof["n_cells"], axis=0).round()
        don = C.donor_ages(n).set_index("donor")
        counts = counts[[d for d in counts.columns if d in don.index]]
        parent = parent_of(prof, ds)
        ok = prof["purity_cell_class"] >= MIN_PURITY
        for par, clusters in parent[ok].groupby(parent[ok]).groups.items():
            if len(clusters) < 2:
                continue
            sub = counts.loc[list(clusters)]
            tot = sub.sum(axis=0)
            donors = [d for d in sub.columns if tot[d] >= MIN_PARENT_CELLS]
            if len(donors) < MIN_DONORS:
                continue
            sub = sub[donors]
            ages = don.loc[donors, "age_pcw"].to_numpy(float)
            lr = C.clr(sub.T).T                      # clusters x donors
            rho = C.spearman_rows(lr.to_numpy(float), ages)
            p, exact = C.spearman_perm_p(rho, ages)
            frac = sub.div(sub.sum(axis=0), axis=1)
            young, old = ages == ages.min(), ages == ages.max()
            for i, c in enumerate(sub.index):
                rows.append({"dataset": ds, "chemistry": chem, "parent": par, "cluster": c,
                             "n_donors": len(donors), "ages": f"{ages.min():g}-{ages.max():g}",
                             "rho_clr_vs_age": rho[i], "perm_p": p[i], "exact": exact,
                             "share_youngest": float(frac.loc[c, young].mean()),
                             "share_oldest": float(frac.loc[c, old].mean()),
                             "mean_share": float(frac.loc[c].mean()),
                             "cells_in_cluster": int(sub.loc[c].sum()),
                             "purity_donor": float(prof.loc[c, "purity_donor"])})
            fl = frac.T.assign(age_pcw=ages).melt(id_vars="age_pcw", var_name="cluster", value_name="share")
            shares.append(fl.assign(dataset=ds, chemistry=chem, parent=par))
    return pd.DataFrame(rows), pd.concat(shares, ignore_index=True) if shares else pd.DataFrame()


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    agree = {ds: C.cluster_label_identity(ds, CLUSTERING[ds], min_cells=C.MIN_CELLS)[0] for ds in C.DATASETS}
    for ds, a in agree.items():
        C.log(f"  {ds} {CLUSTERING[ds]}: {a:.0%} of labels best-match themselves across chemistries")
        for chem in C.CHEMISTRIES:
            out.used(f"{C.ns(ds, chem)}/09_pseudobulk/cluster_{CLUSTERING[ds]}__pseudobulk_counts.csv")
    per, shares = per_stratum(out)
    out.write(per, "cluster_abundance_per_stratum",
              "Per stratum x parent x cluster: Spearman of the cluster's CLR share (within parent) "
              "with donor age, exact permutation p")
    out.write(shares, "cluster_shares_by_donor", "Per stratum x parent: each cluster's share per donor")
    comb = C.combine_chemistries(per[per.dataset.map(lambda d: agree[d] >= C.SHARED_LABEL_IDENTITY)],
                                 ["dataset", "parent", "cluster"], effect="rho_clr_vs_age",
                                 weight="n_donors", labels=("expands", "shrinks"),
                                 carry=("share_youngest", "share_oldest", "n_donors", "ages", "purity_donor"))
    if not comb.empty:
        for ds in C.DATASETS:
            mk = C.top_markers(C.ns(ds, "v2"), f"cluster_{CLUSTERING[ds]}")
            for chem in C.CHEMISTRIES:
                if (C.EXPORTS / C.ns(ds, chem) / "10_markers" / f"top_markers_cluster_{CLUSTERING[ds]}.csv").exists():
                    out.used(f"{C.ns(ds, chem)}/10_markers/top_markers_cluster_{CLUSTERING[ds]}.csv")
            m = comb.dataset == ds
            comb.loc[m, "top_markers"] = comb.loc[m, "cluster"].map(mk).fillna("")
    out.write(comb, "cluster_abundance_combined",
              "Per parent x cluster: v2 x v3 combined age trend of its share within the parent; tier")
    figures(out, comb, shares)

    # ---- findings -----------------------------------------------------------
    f = []
    for ds in C.DATASETS:
        if agree[ds] < C.SHARED_LABEL_IDENTITY:
            f.append(f"**{ds}: not combined** -- {CLUSTERING[ds]} labels are not shared across "
                     f"chemistries (only {agree[ds]:.0%} best-match themselves).")
            continue
        g = comb[comb.dataset == ds] if not comb.empty else comb
        tested = per[(per.dataset == ds)]
        npar = tested.parent.nunique()
        rep = g[g.tier != ""] if len(g) else g
        label = "" if ds == "cortex" else " (parent = class x region)"
        if rep.empty:
            f.append(f"**{ds}{label}: no sub-type share changes with age in both donor sets** "
                     f"({len(g)} clusters in {npar} parents tested).")
            continue
        f.append(f"**{ds}{label}: sub-types whose share within their class changes with age** "
                 f"({len(rep)} of {len(g)} clusters tiered; share youngest -> oldest donor, v2 / v3): "
                 + "; ".join(f"cluster {r.cluster} in {r.parent} {r.direction}"
                             + (f" [{r.top_markers}]" if r.top_markers else "")
                             + f" ({r.share_youngest_v2:.0%}->{r.share_oldest_v2:.0%} / "
                               f"{r.share_youngest_v3:.0%}->{r.share_oldest_v3:.0%}; {r.tier})"
                             for r in rep.head(MAX_LISTED * 2).itertuples())
                 + ("" if len(rep) <= MAX_LISTED * 2 else f" (+{len(rep) - MAX_LISTED * 2} more)") + ".")
    out.summary(
        TITLE,
        "Within each cell class, which sub-types (clusters) expand or shrink with age, "
        "consistently in both donor sets?",
        [f"Clusterings with labels shared across chemistries: cortex {CLUSTERING['cortex']}, human_dev "
         f"{CLUSTERING['human_dev']} (share of labels whose v2 expression profile best matches the "
         "same label in v3: "
         + ", ".join(f"{d} {a:.0%}" for d, a in agree.items()) + ").",
         f"Cells per cluster x donor from 04's cluster sizes and donor fractions; parent = dominant "
         f"class (cortex) or class x dominant region (human_dev); clusters with class purity >= "
         f"{MIN_PURITY}; donors with >= {MIN_PARENT_CELLS} cells in the parent, >= {MIN_DONORS} donors.",
         "Per donor, CLR of counts across the parent's clusters; Spearman with donor age, exact "
         "permutation p; v2 x v3 signed Stouffer (weights sqrt(donors)), BH per dataset, tiered."],
        f,
        ["Clusters were defined on all ages together; a sub-type that exists only late appears as "
         "an expanding cluster, which is the intended reading, but cluster boundaries are not "
         "independent of age.",
         "human_dev regions were dissected differently per donor; class x region parents reduce, "
         "but do not remove, that confound.",
         "5-15 donors per chemistry: a trend is across that many people."],
        ["B6: Milo neighbourhoods (cortex nhoods) test abundance at finer resolution without "
         "fixed cluster boundaries."])


def figures(out: C.Output, comb: pd.DataFrame, shares: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or comb.empty or shares.empty:
        return
    for ds in C.DATASETS:
        top = comb[(comb.dataset == ds) & (comb.tier != "")].head(MAX_LISTED)
        if top.empty:
            continue
        fig, axes = plt.subplots(1, len(top), figsize=(2.6 * len(top), 2.6), squeeze=False)
        for ax, r in zip(axes.flat, top.itertuples()):
            for chem, mk in zip(C.CHEMISTRIES, ["o", "^"]):
                h = shares[(shares.dataset == ds) & (shares.chemistry == chem) & (shares.parent == r.parent)
                           & (shares.cluster == r.cluster)].sort_values("age_pcw")
                ax.plot(h.age_pcw, h.share, marker=mk, lw=0.8, label=chem)
            ax.set_title(f"{r.cluster} in {r.parent}"[:40], fontsize=7)
            ax.set_xlabel("age (pcw)", fontsize=7)
            ax.tick_params(labelsize=6)
        axes.flat[0].set_ylabel("share within parent", fontsize=7)
        axes.flat[0].legend(fontsize=6)
        fig.tight_layout()
        out.figure(fig, f"cluster_abundance_{ds}", f"Sub-types whose within-class share changes with age, {ds}")
        plt.close(fig)


if __name__ == "__main__":
    main()
