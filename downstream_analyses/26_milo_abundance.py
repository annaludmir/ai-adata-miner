#!/usr/bin/env python3
"""26 - Differential abundance over age in Milo neighbourhoods (cortex).

Question: 12 tested sub-type abundance with fixed clusters. Milo's
neighbourhoods -- small, overlapping groups of similar cells stored in the
cortex file -- test it without cluster boundaries and at finer resolution.
Which neighbourhoods grow or shrink with age within their cell class, in both
donor sets, and where do they sit (class, cluster, cycling share)?

Method
  Cells per neighbourhood x donor from stage-1 script 21. Neighbourhood labels
  come from the whole file, so the same neighbourhood exists in both
  chemistries. Abundance = the neighbourhood's cells as a share of the donor's
  cells of the neighbourhood's dominant class (so a class-level shift, 02's
  result, does not count again): log((count + 0.5) / (class cells + 1)).
  Neighbourhoods with >= 20 cells in the stratum and donors with >= 50 cells
  of the class; >= 5 donors. Spearman with donor age, exact permutation p;
  v2 x v3 signed Stouffer (weights sqrt(donors)), BH, tiered. Overlapping
  neighbourhoods are not independent, so counts of tiered neighbourhoods are
  summaries, not independent discoveries (Milo's spatial FDR is not applied).

Inputs (csv_exports/cortex__<chem>/21_milo/):
  nhood_counts_by_donor.csv, nhood_annotation.csv, donor_totals.csv;
  05_confounds/crosstab_age_x_donor.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "26_milo_abundance"
TITLE = "Neighbourhood (Milo) abundance over age within cell classes (cortex)"
MIN_NHOOD_CELLS = 20
MIN_CLASS_CELLS = 50
MIN_DONORS = 5
MAX_LISTED = 8
SHARED_CLUSTERINGS = ["leiden_scVI", "louvain", "cluster_id"]   # labels mean the same cells in v2 and v3


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    rows, ann_all = [], {}
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        base = C.EXPORTS / n / "21_milo"
        if not (base / "nhood_counts_by_donor.csv").exists():
            continue
        out.used(f"{n}/21_milo/nhood_counts_by_donor.csv", f"{n}/21_milo/nhood_annotation.csv",
                 f"{n}/21_milo/donor_totals.csv", f"{n}/05_confounds/crosstab_age_x_donor.csv")
        cnt = pd.read_csv(base / "nhood_counts_by_donor.csv").set_index("nhood")
        ann = pd.read_csv(base / "nhood_annotation.csv").set_index("nhood")
        tot = pd.read_csv(base / "donor_totals.csv").set_index("donor")
        don = C.donor_ages(n).set_index("donor")
        ann_all[(ds, chem)] = ann
        cnt = cnt[[d for d in cnt.columns if d in don.index and d in tot.index]]
        for cls, idx in ann[ann.n_cells >= MIN_NHOOD_CELLS].groupby("dominant_cell_class").groups.items():
            if cls not in tot.columns:
                continue
            donors = [d for d in cnt.columns if tot.loc[d, cls] >= MIN_CLASS_CELLS]
            if len(donors) < MIN_DONORS:
                continue
            ages = don.loc[donors, "age_pcw"].to_numpy(float)
            c = cnt.loc[list(idx), donors].to_numpy(float)
            ab = np.log((c + 0.5) / (tot.loc[donors, cls].to_numpy(float)[None, :] + 1))
            rho = C.spearman_rows(ab, ages)
            p, _ = C.spearman_perm_p(rho, ages)
            share = c / tot.loc[donors, cls].to_numpy(float)[None, :]
            young, old = ages == ages.min(), ages == ages.max()
            for i, nh in enumerate(idx):
                rows.append({"dataset": ds, "chemistry": chem, "nhood": nh, "cell_class": cls,
                             "n_donors": len(donors), "rho_vs_age": rho[i], "perm_p": p[i],
                             "share_youngest": float(share[i, young].mean()), "share_oldest": float(share[i, old].mean())})
    if not rows:
        out.summary(TITLE, "Skipped: no Milo neighbourhood counts.",
                    ["Needs 21_milo/ from stage-1 script 21 (cortex only: human_dev stores no neighbourhoods)."],
                    ["**Not run**: script 21 is new; re-run stage 1 (slurm_01_metadata.sh), then this analysis."], [], [])
        return
    per = pd.DataFrame(rows)
    comb = C.combine_chemistries(per, ["dataset", "nhood"], effect="rho_vs_age", weight="n_donors",
                                 labels=("expands with age", "shrinks with age"),
                                 carry=("rho_vs_age", "share_youngest", "share_oldest", "cell_class"))
    if not comb.empty:
        ann = ann_all.get((comb.dataset.iloc[0], "v2"))
        extra = [c for c in ann.columns if c.startswith("dominant_") and c != "dominant_cell_class"] + \
            [c for c in ("mean_age_pcw", "cycling_share", "dominant_class_share") if c in ann.columns]
        comb = comb.join(ann[extra], on="nhood")
    summ = pd.DataFrame()
    if not comb.empty:
        comb["cell_class"] = comb["cell_class_v2"].fillna(comb["cell_class_v3"])
        summ = (comb.groupby("cell_class")
                .apply(lambda g: pd.Series({"n_nhoods": len(g), "expanding": int(((g.tier != "") & (g.stouffer_z > 0)).sum()),
                                            "shrinking": int(((g.tier != "") & (g.stouffer_z < 0)).sum())}),
                       include_groups=False).reset_index())
    out.write(per, "nhood_age_trends_per_stratum", "Per stratum x neighbourhood: within-class abundance vs donor age")
    out.write(comb, "nhood_age_trends_combined", "Per neighbourhood: v2 x v3 combined; tier; annotation")
    out.write(summ, "nhood_summary_by_class", "Per class: neighbourhoods tested, tiered expanding and shrinking")

    f = []
    if len(summ):
        f.append("**Neighbourhoods changing with age within their class** (tested / expanding / shrinking, "
                 "tiered): " + "; ".join(f"{r.cell_class} {r.n_nhoods} / {r.expanding} / {r.shrinking}"
                                         for r in summ.itertuples()) + ".")
        # a clustering shared across chemistries labels neighbourhoods consistently
        pref = [f"dominant_{c}" for c in SHARED_CLUSTERINGS]
        cl = [c for c in pref if c in comb.columns]
        rep = comb[comb.tier != ""]
        if cl and len(rep):
            by = rep.groupby([cl[0], "direction"]).size().unstack(fill_value=0)
            f.append(f"**By {cl[0].replace('dominant_', '')} cluster** (tiered expanding / shrinking neighbourhoods): "
                     + "; ".join(f"{k}: {int(r.get('expands with age', 0))} / {int(r.get('shrinks with age', 0))}"
                                 for k, r in by.iterrows()) + ".")
        if "cycling_share" in comb:
            rep = comb[comb.tier != ""]
            if len(rep):
                f.append("**Cycling share of changing neighbourhoods** (median): expanding "
                         f"{rep[rep.stouffer_z > 0].cycling_share.median():.0%}, shrinking "
                         f"{rep[rep.stouffer_z < 0].cycling_share.median():.0%}, all tested {comb.cycling_share.median():.0%}.")
    out.summary(
        TITLE,
        "Which Milo neighbourhoods grow or shrink with age within their cell class, in both donor sets, and "
        "where do they sit?",
        [f"Cells per neighbourhood x donor (stage-1 script 21); within-class share log((count + 0.5) / (class "
         f"cells + 1)); neighbourhoods >= {MIN_NHOOD_CELLS} cells, donors >= {MIN_CLASS_CELLS} class cells, >= "
         f"{MIN_DONORS} donors; Spearman with age, exact permutation; v2 x v3 signed Stouffer, BH, tiered."],
        f,
        ["Neighbourhoods overlap, so tests are correlated and counts of significant neighbourhoods overstate "
         "independent findings; Milo's spatial FDR is not applied.",
         "The neighbourhood graph was built on all cells together, so neighbourhoods are shared across chemistries "
         "only as far as the graph mixed them; v2-only or v3-only neighbourhoods are tested in one chemistry."],
        ["Run Milo's negative-binomial GLM with spatial FDR from the stored graph for a full analysis."])


if __name__ == "__main__":
    main()
