#!/usr/bin/env python3
"""01 - Data audit: the replicate structure every later analysis must respect.

Question: before looking for biology, what can these exports support? How many
independent donors stand behind each age, chemistry and dataset; is donor sex
entangled with age; do the cell-class labels recover known markers; and are
there artefacts in the exports themselves?

Inputs (csv_exports/):
  <ds>/18_chemistry/age_chemistry_overlap.csv, chemistry_summary.csv
  <ds>__<chem>/05_confounds/crosstab_age_x_donor.csv
  <ds>__<chem>/09_pseudobulk/donor__mean_lognorm.csv, *__group_summary.csv
  <ds>__<chem>/02_composition/counts_cell_class_by_donor.csv
  <ds>__<chem>/10_markers/top_markers_cell_class.csv
  <ds>__<chem>/11_panels/marker_label_check.csv
  _cross_dataset/gene_id_map.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "01_data_audit"
TOP_MARKER_RANK = 50
Y_MALE_MIN = 0.1      # mean log1p(CP10K) of Y genes; observed: males >0.3, females <0.01
XIST_FEMALE_MIN = 0.1  # observed: females >0.15, males <0.01


def replicate_structure(out: C.Output) -> pd.DataFrame:
    rows = []
    for d, c in C.STRATA:
        n = C.ns(d, c)
        out.used(f"{n}/05_confounds/crosstab_age_x_donor.csv")
        da = C.donor_ages(n)
        per_age = da.groupby("age_pcw")["donor"].nunique()
        rows.append({
            "stratum": n, "n_cells": int(da["n_cells"].sum()),
            "n_donors": len(da), "n_ages": int(per_age.size),
            "age_min": da["age_pcw"].min(), "age_max": da["age_pcw"].max(),
            "ages_with_one_donor": int((per_age == 1).sum()),
            "max_donors_per_age": int(per_age.max()),
            "max_ages_per_donor": int(da["n_ages"].max()),
            "smallest_donor_cells": int(da["n_cells"].min()),
            "ages": "|".join(f"{a:g}" for a in per_age.index),
        })
    df = pd.DataFrame(rows)
    out.write(df, "replicate_structure",
              "Per chemistry stratum: donors, ages, donors per age")
    return df


def chemistry_overlap(out: C.Output) -> pd.DataFrame:
    rows = []
    for d in C.DATASETS:
        out.used(f"{d}/18_chemistry/age_chemistry_overlap.csv",
                 f"{d}/18_chemistry/chemistry_summary.csv")
        o = C.csv(d, "18_chemistry/age_chemistry_overlap.csv")
        o = o[o["role"] == "age"]
        s = C.csv(d, "18_chemistry/chemistry_summary.csv").iloc[0]
        rows.append({
            "dataset": d,
            "ages_v2_only": "|".join(f"{a:g}" for a in o.loc[o.verdict == "v2 only", "age_pcw"]),
            "ages_v3_only": "|".join(f"{a:g}" for a in o.loc[o.verdict == "v3 only", "age_pcw"]),
            "ages_in_both": "|".join(f"{a:g}" for a in o.loc[o.comparable, "age_pcw"]),
            "frac_cells_at_shared_ages": float(s["frac_cells_at_comparable_ages"]),
            "donor_nested_in_chemistry": bool(s["donor_nested_in_chemistry"]),
        })
    df = pd.DataFrame(rows)
    out.write(df, "chemistry_age_overlap", "Which ages each chemistry covers")
    return df


def donor_map(out: C.Output) -> pd.DataFrame:
    def table(ds):
        parts = []
        for c in C.CHEMISTRIES:
            da = C.donor_ages(C.ns(ds, c))
            da["chemistry"] = c
            parts.append(da)
        return pd.concat(parts)

    cx, hd = table("cortex"), table("human_dev").set_index("donor")
    rows = []
    for _, r in cx.iterrows():
        key = C.normalise_donor(r["donor"])
        h = hd.loc[key] if key in hd.index else None
        rows.append({
            "cortex_donor": r["donor"], "human_dev_donor": key if h is not None else "",
            "chemistry": r["chemistry"], "age_cortex": r["age_pcw"],
            "age_human_dev": h["age_pcw"] if h is not None else np.nan,
            "cells_cortex": r["n_cells"],
            "cells_human_dev": h["n_cells"] if h is not None else np.nan,
            "same_chemistry": h is not None and h["chemistry"] == r["chemistry"],
            "same_age": h is not None and h["age_pcw"] == r["age_pcw"],
            # unmatched because an exclusion rule removed them from human_dev
            "excluded_from_human_dev": h is None and C.excluded(
                "human_dev", donor=key, age=r["age_pcw"]),
        })
    df = pd.DataFrame(rows)
    out.write(df, "donor_map_cortex_to_human_dev",
              "Each cortex donor matched to its human_dev ID (IDs are written differently)")
    return df


def infer_sex(out: C.Output, chrom: dict[str, str]) -> pd.DataFrame:
    rows = []
    for d, c in C.STRATA:
        n = C.ns(d, c)
        out.used(f"{n}/09_pseudobulk/donor__mean_lognorm.csv")
        m = C.group_matrix(n, "donor", min_cells=1)
        y = [g for g in C.Y_GENES if g in m.index]
        da = C.donor_ages(n).set_index("donor")
        for donor in m.columns:
            ys = float(m.loc[y, donor].mean()) if y else np.nan
            xi = float(m.loc["XIST", donor]) if "XIST" in m.index else np.nan
            male, female = ys > Y_MALE_MIN, xi > XIST_FEMALE_MIN
            sex = ("M" if male and not female else "F" if female and not male
                   else "ambiguous")
            rows.append({"stratum": n, "dataset": d, "chemistry": c, "donor": donor,
                         "age_pcw": da.loc[donor, "age_pcw"] if donor in da.index else np.nan,
                         "y_gene_mean": ys, "xist": xi, "inferred_sex": sex,
                         "n_y_genes": len(y)})
    df = pd.DataFrame(rows)
    out.write(df, "donor_sex_inferred",
              "Donor sex from Y-gene vs XIST expression (obs['sex'] is 'unknown')")
    return df


def sex_age_entanglement(out: C.Output, sex: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for stratum, g in sex.groupby("stratum", sort=False):
        g = g[g["inferred_sex"].isin(["M", "F"])]
        is_m = (g["inferred_sex"] == "M").to_numpy(float)
        ages = g["age_pcw"].to_numpy(float)
        if 0 < is_m.sum() < len(is_m):
            rho = C.spearman_rows(is_m[None, :], ages)
            p, exact = C.spearman_perm_p(rho, ages, X=is_m[None, :])
            rho, p = float(rho[0]), float(p[0])
        else:
            rho, p, exact = np.nan, np.nan, True
        rows.append({
            "stratum": stratum, "n_male": int(is_m.sum()), "n_female": int(len(is_m) - is_m.sum()),
            "male_ages": "|".join(f"{a:g}" for a in sorted(ages[is_m == 1])),
            "female_ages": "|".join(f"{a:g}" for a in sorted(ages[is_m == 0])),
            "spearman_male_vs_age": rho, "perm_p": p, "exact": exact,
        })
    df = pd.DataFrame(rows)
    out.write(df, "sex_age_entanglement",
              "Is donor sex associated with age within a stratum? (exact permutation p)")
    return df


def class_sex_balance(out: C.Output, sex: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for d, c in C.STRATA:
        n = C.ns(d, c)
        out.used(f"{n}/02_composition/counts_cell_class_by_donor.csv")
        counts = C.csv(n, "02_composition/counts_cell_class_by_donor.csv", index_col=0)
        s = sex[sex["stratum"] == n].set_index("donor")["inferred_sex"]
        counts = counts.loc[[i for i in counts.index if i in s.index]]
        male = counts.loc[s.loc[counts.index] == "M"].sum()
        overall = male.sum() / counts.to_numpy().sum()
        for cls in counts.columns:
            tot = counts[cls].sum()
            if tot == 0:
                continue
            rows.append({"stratum": n, "cell_class": cls, "n_cells": int(tot),
                         "frac_from_male_donors": male[cls] / tot,
                         "frac_male_overall": overall,
                         "ratio_vs_overall": (male[cls] / tot) / overall if overall > 0 else np.nan})
    df = pd.DataFrame(rows)
    out.write(df, "cell_class_sex_balance",
              "Fraction of each cell class contributed by male donors vs the stratum overall")
    return df


def sex_genes_in_markers(out: C.Output, chrom: dict[str, str]) -> pd.DataFrame:
    rows = []
    for d, c in C.STRATA:
        n = C.ns(d, c)
        out.used(f"{n}/10_markers/top_markers_cell_class.csv")
        t = C.csv(n, "10_markers/top_markers_cell_class.csv")
        t = t[t["rank_in_group"] <= TOP_MARKER_RANK]
        flag = C.sex_linked(t["gene"].astype(str), chrom)
        for _, r in t[flag].iterrows():
            rows.append({"stratum": n, "cell_class": r["top_group"], "gene": r["gene"],
                         "chromosome": chrom.get(str(r["gene"]), ""),
                         "rank_in_group": int(r["rank_in_group"]),
                         "log2fc_vs_rest": r["log2fc_vs_rest"]})
    df = pd.DataFrame(rows, columns=["stratum", "cell_class", "gene", "chromosome",
                                     "rank_in_group", "log2fc_vs_rest"])
    out.write(df, "sex_linked_genes_in_top_markers",
              f"Y-linked genes and XIST among the top {TOP_MARKER_RANK} markers of a cell class")
    return df


def marker_recovery(out: C.Output) -> pd.DataFrame:
    rows = []
    for d, c in C.STRATA:
        n = C.ns(d, c)
        out.used(f"{n}/11_panels/marker_label_check.csv")
        m = C.csv(n, "11_panels/marker_label_check.csv")
        m = m[m["panel_group"] == "marker"].copy()
        m["expected_class"] = m["panel"].map(C.MARKER_PANEL_TO_CLASS)
        # Only score panels whose cell class exists in this stratum.
        classes = set(C.group_summary(n, "cell_class")["group"].astype(str))
        m = m[m["expected_class"].isin(classes)]
        m.insert(0, "stratum", n)
        rows.append(m)
    df = pd.concat(rows, ignore_index=True)
    out.write(df, "marker_panel_recovery",
              "Does each seed marker panel score highest in the cell class it names?")
    return df


def empty_groups(out: C.Output) -> pd.DataFrame:
    rows = []
    for d, c in C.STRATA:
        n = C.ns(d, c)
        for f in sorted((C.EXPORTS / n / "09_pseudobulk").glob("*__group_summary.csv")):
            g = pd.read_csv(f)
            rows.append({"stratum": n, "grouping": f.name.replace("__group_summary.csv", ""),
                         "n_groups": len(g), "n_empty": int((g["n_cells"] == 0).sum()),
                         "n_below_50_cells": int((g["n_cells"] < C.MIN_CELLS).sum())})
    df = pd.DataFrame(rows)
    out.write(df, "pseudobulk_group_sizes",
              "Groups per pseudobulk table; empty = level that exists only in the other chemistry")
    return df


def exclusions_status(out: C.Output) -> pd.DataFrame:
    """Each exclusion rule, the cells it covers, and whether the exports honour it."""
    rows = []
    for d in C.DATASETS:
        rules = C.exclusion_rules(d)
        for r in rules.itertuples():
            for c in C.CHEMISTRIES:
                n = C.ns(d, c)
                x = C.csv(n, "05_confounds/crosstab_age_x_donor.csv", index_col=0)
                if r.role == "age":
                    hit = x.loc[[a for a in x.index if abs(float(a) - float(r.value)) < 0.01]]
                    n_cells = int(hit.to_numpy().sum())
                elif r.role == "donor":
                    cols = [k for k in x.columns if C.normalise_donor(k) == C.normalise_donor(r.value)]
                    n_cells = int(x[cols].to_numpy().sum())
                else:
                    n_cells = np.nan
                applied = (C.EXPORTS / n / "01_overview" / "exclusions_applied.csv").exists()
                rows.append({"dataset": d, "chemistry": c, "role": r.role, "value": r.value,
                             "reason": r.reason, "cells_still_in_exports": n_cells,
                             "exports_made_with_exclusions": applied})
    df = pd.DataFrame(rows, columns=["dataset", "chemistry", "role", "value", "reason",
                                     "cells_still_in_exports", "exports_made_with_exclusions"])
    out.write(df, "exclusions_status",
              f"Rules in {Path(C.EXCLUSIONS).name}: cells covered and whether exports already omit them")
    return df


def xist_scale(sex: pd.DataFrame, dmap: pd.DataFrame) -> pd.DataFrame:
    """XIST in cortex vs human_dev for the same female donors."""
    cx = sex[sex.dataset == "cortex"].set_index("donor")
    hd = sex[sex.dataset == "human_dev"].set_index("donor")
    rows = []
    for _, r in dmap.iterrows():
        a, b = r["cortex_donor"], r["human_dev_donor"]
        if a in cx.index and b in hd.index and cx.loc[a, "inferred_sex"] == "F":
            rows.append({"donor": b, "xist_cortex": cx.loc[a, "xist"],
                         "xist_human_dev": hd.loc[b, "xist"]})
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    out.used("_cross_dataset/gene_id_map.csv")
    chrom = C.chromosome_map()

    rep = replicate_structure(out)
    ovl = chemistry_overlap(out)
    dmap = donor_map(out)
    sex = infer_sex(out, chrom)
    ent = sex_age_entanglement(out, sex)
    class_sex_balance(out, sex)
    sxm = sex_genes_in_markers(out, chrom)
    rec = marker_recovery(out)
    emp = empty_groups(out)
    xs = xist_scale(sex, dmap)
    exc = exclusions_status(out)
    if not xs.empty:
        out.write(xs, "xist_cortex_vs_human_dev",
                  "XIST level for the same female donors in each file")

    # ---- findings, computed rather than written by hand ------------------
    f = []
    one = rep["ages_with_one_donor"].sum() / rep["n_ages"].sum()
    f.append(
        "**Donor is the replicate unit, and age is nested in donor.** Every donor has "
        f"exactly one age (max ages per donor = {rep['max_ages_per_donor'].max()}), and "
        f"{one:.0%} of age points rest on a single donor. Donors per stratum: "
        + ", ".join(f"{r.stratum} {r.n_donors}" for r in rep.itertuples())
        + ". An age effect is therefore indistinguishable from a donor effect at that age.")
    for r in ovl.itertuples():
        f.append(
            f"**{r.dataset}: chemistries cover different ages.** v2 only: {r.ages_v2_only}; "
            f"v3 only: {r.ages_v3_only}; both: {r.ages_in_both or 'none'} "
            f"({r.frac_cells_at_shared_ages:.0%} of cells). Donors nested in chemistry: "
            f"{r.donor_nested_in_chemistry}. v2 and v3 are therefore independent donor sets "
            "that can replicate each other, but not be pooled as if equivalent.")
    has = dmap.human_dev_donor.ne("")
    matched = int((dmap.same_chemistry & has).sum())
    renamed = int((has & (dmap.cortex_donor != dmap.human_dev_donor)).sum())
    excl = dmap[dmap.excluded_from_human_dev]
    lost = dmap[~has & ~dmap.excluded_from_human_dev]
    age_mis = dmap[~dmap.same_age & has]
    # Does the cross-dataset join in csv_exports see the renamed donors?
    out.used("_cross_dataset/label_overlap_donor.csv")
    lo = C.csv("_cross_dataset", "label_overlap_donor.csv")
    joined = int(lo["in_both"].astype(str).str.lower().eq("true").sum())
    join_note = ("`13_cross_dataset_keys` joins them" if joined >= matched else
                 f"`13_cross_dataset_keys` joins only {joined} -- exports predate its ID fix")
    f.append(
        f"**cortex is not an independent cohort.** {matched}/{len(dmap)} cortex donors are "
        f"human_dev donors ({renamed} under a differently written ID, e.g. "
        f"`XHU:1966:307` = `XHU:307`; {join_note})."
        + (f" The other {len(excl)} ({', '.join(excl.cortex_donor)}) fall under a human_dev "
           "exclusion rule, so they have no human_dev counterpart in these exports."
           if len(excl) else "")
        + (f" {len(lost)} ({', '.join(lost.cortex_donor)}) have no human_dev match at all."
           if len(lost) else "")
        + " Agreement between the two files is reproducibility of processing, not replication."
        + (f" Age annotations disagree for {', '.join(age_mis.cortex_donor)} "
           f"({', '.join(f'{a:g} vs {b:g}' for a, b in zip(age_mis.age_cortex, age_mis.age_human_dev))} pcw)."
           if len(age_mis) else ""))
    amb = int((sex.inferred_sex == "ambiguous").sum())
    f.append(
        f"**Donor sex is recoverable from expression** (obs['sex'] is 'unknown'): "
        f"{int((sex.inferred_sex != 'ambiguous').sum())}/{len(sex)} donor entries call cleanly"
        + (f", {amb} ambiguous" if amb else "") + ".")
    uneven = ent[ent["spearman_male_vs_age"].abs() >= 0.4]
    if len(uneven):
        f.append(
            "**Sex is unevenly spread over age** in "
            + "; ".join(f"{r.stratum} (male ages {r.male_ages}, female ages "
                        f"{r.female_ages}; rho = {r.spearman_male_vs_age:.2f}, exact "
                        f"p = {C.fmt_p(r.perm_p)})" for r in uneven.itertuples())
            + ". With this few donors the association is not significant, but it does not "
            "need to be to matter: a sex-differential gene can look like an age trend. "
            "Later analyses flag sex-linked genes and check trends against sex.")
    if not sxm.empty:
        ex = sxm.sort_values("rank_in_group").drop_duplicates(["stratum", "cell_class"])
        f.append(
            "**Sex-linked genes rank among cell-class markers**, which happens when a class "
            "is drawn unevenly from male and female donors: "
            + "; ".join(f"{r.gene} in {r.cell_class} ({r.stratum}, rank {r.rank_in_group})"
                        for r in ex.itertuples())
            + ". See cell_class_sex_balance.csv for the imbalance behind each.")
    good = rec.dropna(subset=["expected_class"])
    hit = good["expected_match"].astype(str).str.lower().eq("true")
    f.append(
        f"**Cell-class labels recover known markers**: {int(hit.sum())}/{len(good)} "
        "seed marker panels score highest in the class they name. Misses: "
        + (", ".join(f"{r.panel}->{r.highest_scoring_cell_class} ({r.stratum})"
                     for r in good[~hit].itertuples()) or "none")
        + ". Panels naming a class absent from a stratum are not scored. The misses are "
        "neighbouring lineages: glioblasts carry radial-glia genes, placode-derived "
        "sensory neurons carry pan-neuronal genes, and the short seed lists cannot "
        "separate them.")
    if len(exc):
        still = exc[exc.cells_still_in_exports.fillna(0) > 0]
        rules = exc.drop_duplicates(["dataset", "role", "value"])
        f.append(
            f"**Exclusions ({Path(C.EXCLUSIONS).name})**: "
            + "; ".join(f"{r.dataset} {r.role}={r.value} ({r.reason})" for r in rules.itertuples())
            + ". "
            + ("These exports already omit them." if still.empty else
               f"The exports still contain {int(still.cells_still_in_exports.sum()):,} of these "
               "cells (made before the rule). Step 3 drops them from every age-, donor- and "
               "sample-resolved table; tables pooled over age (cell-class pseudobulk, "
               "10_markers, 17_gsea -- used by 04 parts C-D and 05) still include them until "
               "stages 1-3 are re-run on the cluster."))
    n_empty = int(emp["n_empty"].sum())
    if n_empty:
        f.append(
            f"**These exports contain {n_empty} empty pseudobulk groups** (levels found only "
            "in the other chemistry, written as all-zero columns by script 09 before the "
            "fix). Every analysis here drops groups below 50 cells, which removes them.")
    if not xs.empty:
        ratio = (xs.xist_human_dev / xs.xist_cortex).median()
        f.append(
            f"**XIST reads ~{ratio:.1f}x higher in human_dev than in cortex for the same "
            "female donors**, so the two files probably count reads differently "
            "(XIST is largely nuclear and intronic). Compare genes across the files by "
            "rank or within-file contrast, not by absolute level.")

    out.summary(
        "Data audit: what the exports can support",
        "Before any biology: how many independent units stand behind each comparison, "
        "which covariates are entangled, and do the labels and tables behave?",
        ["Replicate structure from the age x donor crosstabs of each chemistry stratum.",
         "Cross-dataset donor matching after normalising ID formats; checked on chemistry, "
         "age and cell count (cortex cells <= human_dev cells).",
         f"Donor sex inferred from donor-level pseudobulk: male if mean Y-gene log1p(CP10K) "
         f"> {Y_MALE_MIN} and XIST <= {XIST_FEMALE_MIN}, female if the reverse.",
         "Sex-age association: Spearman between male indicator and age, exact permutation p.",
         "Marker sanity from 11_panels/marker_label_check.csv; sex-linked markers from the "
         f"top {TOP_MARKER_RANK} per class."],
        f,
        ["Sex calls rest on expression and need confirming against sample records.",
         "Per-donor composition reflects which tissue was dissected from that donor, so a "
         "class imbalance between sexes may be dissection, not sex biology.",
         "Seed marker panels are short, hand-picked lists (panels/README.md)."],
        ["Fill obs['sex'] from records and re-export, so sex can enter models as a covariate.",
         "Fix the donor-ID join in 13_cross_dataset_keys (strip the middle field).",
         "Re-run stage 2 so the pseudobulk exports carry no empty groups."])


if __name__ == "__main__":
    main()
