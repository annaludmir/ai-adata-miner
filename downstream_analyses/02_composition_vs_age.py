#!/usr/bin/env python3
"""02 - Does cell-class composition change with age, replicated across donor sets?

Question: as development proceeds, which cell classes expand or shrink -- and
does the trend hold in both chemistries, which are independent donor sets?

Method: cells are summed per donor (the replicate unit) within one tissue scope,
turned into a centred log-ratio (composition-aware: fractions cannot move
independently), and correlated with age by Spearman with an exact permutation
p. The v2 and v3 trends are then combined by a signed Stouffer test, so a trend
that flips direction between donor sets cancels rather than passes.

Scopes: cortex uses all cells (the file is the EMX1 dorsal lineage). human_dev
is analysed within one region at a time, because which region was dissected
from a donor changes with age -- whole-brain composition would mostly measure
dissection, not development.

Inputs (csv_exports/):
  <ds>__<chem>/02_composition/counts_cell_class_by_sample.csv
  <ds>__<chem>/01_overview/sample_summary.csv
  results/01_data_audit/donor_sex_inferred.csv  (run 01 first)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "02_composition_vs_age"
MIN_DONORS = 5          # per chemistry, within a scope
MIN_CLASS_FRACTION = 0.005  # rarer classes are merged into "Other" before CLR
SEX_RHO = 0.5           # |Spearman| with the male indicator counted as tracking sex


def donor_counts(n: str, region: str | None) -> pd.DataFrame:
    counts = C.csv(n, "02_composition/counts_cell_class_by_sample.csv", index_col=0)
    meta = C.csv(n, "01_overview/sample_summary.csv").set_index("sample")
    meta = meta.loc[[s for s in counts.index if s in meta.index]]
    if region is not None:
        meta = meta[meta["region"] == region]
    ds = n.split("__")[0]
    keep = np.array([not C.excluded(ds, sample=s, donor=r.donor, age=r.age_pcw)
                     for s, r in meta.iterrows()], dtype=bool)
    meta = meta.loc[keep]
    counts = counts.loc[meta.index]
    by_donor = counts.groupby(meta["donor"]).sum()
    by_donor = by_donor.loc[by_donor.sum(axis=1) >= C.MIN_DONOR_CELLS]
    return by_donor.loc[:, by_donor.sum() > 0]


def scopes() -> list[tuple[str, str | None]]:
    out = [("cortex", None)]
    # Count donors the way the analysis will see them: after exclusions and
    # the per-donor cell minimum, in each chemistry.
    regions = set()
    for c in C.CHEMISTRIES:
        regions |= set(C.csv(C.ns("human_dev", c), "01_overview/sample_summary.csv")["region"])
    shared = sorted(r for r in regions
                    if all(len(donor_counts(C.ns("human_dev", c), r)) >= MIN_DONORS
                           for c in C.CHEMISTRIES))
    return out + [("human_dev", r) for r in shared]


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    sex_path = C.RESULTS / "01_data_audit" / "donor_sex_inferred.csv"
    sex = pd.read_csv(C.require(sex_path))
    out.used("results/01_data_audit/donor_sex_inferred.csv")

    trend_rows, frac_rows = [], []
    for ds, region in scopes():
        scope = region or "all cells"
        per_chem = {}
        for c in C.CHEMISTRIES:
            n = C.ns(ds, c)
            out.used(f"{n}/02_composition/counts_cell_class_by_sample.csv",
                     f"{n}/01_overview/sample_summary.csv")
            per_chem[c] = donor_counts(n, region)
        # One class set per scope so v2 and v3 CLRs are on the same parts.
        pooled = sum(per_chem[c].sum() for c in C.CHEMISTRIES)
        frac = pooled / pooled.sum()
        keep = sorted(frac[frac >= MIN_CLASS_FRACTION].index)
        rare = [k for k in frac.index if k not in keep]

        for c in C.CHEMISTRIES:
            n = C.ns(ds, c)
            cnt = per_chem[c].reindex(columns=frac.index, fill_value=0)
            parts = cnt[keep].copy()
            if rare:
                parts["Other"] = cnt[rare].sum(axis=1)
            ages = C.donor_ages(n).set_index("donor")["age_pcw"]
            sx = sex[sex.stratum == n].set_index("donor")["inferred_sex"]
            parts = parts.loc[[d for d in parts.index if d in ages.index]]
            if len(parts) < MIN_DONORS:
                C.log(f"  {ds} {c} {scope}: {len(parts)} donors (<{MIN_DONORS}), skipped")
                continue
            a = ages.loc[parts.index].to_numpy(float)
            male = (sx.reindex(parts.index) == "M").to_numpy(float)
            z = C.clr(parts)
            fr = parts.div(parts.sum(axis=1), axis=0)
            for donor in parts.index:
                for cls in parts.columns:
                    frac_rows.append({"dataset": ds, "scope": scope, "chemistry": c,
                                      "donor": donor, "age_pcw": ages[donor],
                                      "sex": sx.get(donor, ""), "cell_class": cls,
                                      "n_cells": int(parts.loc[donor, cls]),
                                      "fraction": fr.loc[donor, cls],
                                      "clr": z.loc[donor, cls]})
            X = z.T.to_numpy()
            rho = C.spearman_rows(X, a)
            p, exact = C.spearman_perm_p(rho, a, X=X)
            has_both = 0 < male.sum() < len(male)
            rho_sex = C.spearman_rows(X, male) if has_both else np.full(len(rho), np.nan)
            for i, cls in enumerate(z.columns):
                trend_rows.append({
                    "dataset": ds, "scope": scope, "chemistry": c, "cell_class": cls,
                    "n_donors": len(parts), "age_range": f"{a.min():g}-{a.max():g}",
                    "spearman_rho_age": rho[i], "perm_p": p[i], "exact": exact,
                    "spearman_rho_male": rho_sex[i],
                    "mean_fraction": float(fr[cls].mean())})

    trends = pd.DataFrame(trend_rows)
    fracs = pd.DataFrame(frac_rows)
    out.write(fracs, "donor_composition", "Per donor: cell-class count, fraction and CLR")
    out.write(trends, "composition_age_trends",
              "Per stratum: Spearman of CLR vs age across donors, exact permutation p")

    # ---- replication: combine v2 and v3 per (dataset, scope, class) --------
    rep_rows = []
    for (ds, scope, cls), g in trends.groupby(["dataset", "scope", "cell_class"], sort=False):
        g = g.set_index("chemistry")
        if not set(C.CHEMISTRIES) <= set(g.index):
            continue
        r2, r3 = g.loc["v2"], g.loc["v3"]
        Z, pc = C.signed_stouffer([np.array([r2.spearman_rho_age]), np.array([r3.spearman_rho_age])],
                                  [np.array([r2.perm_p]), np.array([r3.perm_p])],
                                  [np.sqrt(r2.n_donors), np.sqrt(r3.n_donors)])
        rep_rows.append({
            "dataset": ds, "scope": scope, "cell_class": cls,
            "rho_v2": r2.spearman_rho_age, "p_v2": r2.perm_p, "n_v2": int(r2.n_donors),
            "rho_v3": r3.spearman_rho_age, "p_v3": r3.perm_p, "n_v3": int(r3.n_donors),
            "same_direction": np.sign(r2.spearman_rho_age) == np.sign(r3.spearman_rho_age),
            "stouffer_z": float(Z[0]), "combined_p": float(pc[0]),
            "rho_male_v2": r2.spearman_rho_male, "rho_male_v3": r3.spearman_rho_male,
        })
    rep = pd.DataFrame(rep_rows)
    rep["combined_q"] = np.nan
    for key, idx in rep.groupby(["dataset", "scope"]).groups.items():
        rep.loc[idx, "combined_q"] = C.bh(rep.loc[idx, "combined_p"])
    rep["tier"] = C.replication_tier(rep.rho_v2, rep.p_v2, rep.rho_v3, rep.p_v3,
                                     rep.combined_q)
    rep["replicated"] = rep["tier"] == "replicated"
    rep["direction"] = np.where(rep.stouffer_z > 0, "increases with age", "decreases with age")
    rep = rep.sort_values(["dataset", "scope", "combined_p"])
    out.write(rep, "composition_trends_replicated",
              "v2 and v3 trends combined (signed Stouffer); tier: replicated = both "
              "chemistries nominal + combined q < 0.05, supported = combined only")

    figures(out, fracs)

    # ---- findings ----------------------------------------------------------
    f = []
    for (ds, scope), g in rep.groupby(["dataset", "scope"], sort=False):
        hits = g[g.replicated]
        sup = g[g.tier == "supported"]
        sup_txt = ("" if sup.empty else " Supported by the combined test but weak in one "
                   "chemistry: " + "; ".join(
                       f"{r.cell_class} {'up' if r.stouffer_z > 0 else 'down'} "
                       f"(rho v2 {r.rho_v2:+.2f}, v3 {r.rho_v3:+.2f})" for r in sup.itertuples())
                   + ".")
        n2, n3 = int(g.n_v2.iloc[0]), int(g.n_v3.iloc[0])
        if hits.empty:
            best = g.iloc[0]
            f.append(f"**{ds} / {scope}** ({n2} + {n3} donors): no class trend replicates "
                     f"in both chemistries. Strongest: {best.cell_class} (rho v2 "
                     f"{best.rho_v2:+.2f}, v3 {best.rho_v3:+.2f}, combined p = "
                     f"{C.fmt_p(best.combined_p)}).{sup_txt}")
            continue
        parts = []
        for r in hits.itertuples():
            # Sex could only explain a trend replicated in both donor sets if it
            # tracks the class the same way in both; one-sided tracking is the
            # sex/age alignment of that stratum's donors, which the other
            # stratum already controls for.
            sex_note = ""
            if (np.isfinite(r.rho_male_v2) and np.isfinite(r.rho_male_v3)
                    and min(abs(r.rho_male_v2), abs(r.rho_male_v3)) >= SEX_RHO
                    and np.sign(r.rho_male_v2) == np.sign(r.rho_male_v3)):
                sex_note = (f" -- also tracks donor sex in both chemistries "
                            f"(rho {r.rho_male_v2:+.2f}, {r.rho_male_v3:+.2f}), read with care")
            parts.append(f"{r.cell_class} {'up' if r.stouffer_z > 0 else 'down'} "
                         f"(rho v2 {r.rho_v2:+.2f}, v3 {r.rho_v3:+.2f}, "
                         f"q = {C.fmt_p(r.combined_q)}){sex_note}")
        f.append(f"**{ds} / {scope}** ({n2} + {n3} donors), replicated: "
                 + "; ".join(parts) + "." + sup_txt)
    one_sided = rep[rep.replicated & (np.maximum(rep.rho_male_v2.abs(), rep.rho_male_v3.abs()) >= SEX_RHO)
                    & ~((np.minimum(rep.rho_male_v2.abs(), rep.rho_male_v3.abs()) >= SEX_RHO)
                        & (np.sign(rep.rho_male_v2) == np.sign(rep.rho_male_v3)))]
    if len(one_sided):
        f.append("**Sex does not explain the replicated trends.** Where a class tracks donor "
                 "sex, it does so in one chemistry only, or in opposite directions (e.g. "
                 + "; ".join(f"{r.dataset}/{r.scope} {r.cell_class}: rho with male "
                             f"{r.rho_male_v2:+.2f} in v2, {r.rho_male_v3:+.2f} in v3"
                             for r in one_sided.head(3).itertuples())
                 + "), while the age trend keeps its direction in both.")
    flips = rep[~rep.same_direction & (np.minimum(rep.p_v2, rep.p_v3) < 0.05)]
    if len(flips):
        f.append("**Trends that flip between chemistries** (significant in one, opposite sign "
                 "in the other) -- not replicated, possibly donor- or window-specific: "
                 + "; ".join(f"{r.dataset}/{r.scope} {r.cell_class}" for r in flips.itertuples())
                 + ".")

    out.summary(
        "Cell-class composition across age",
        "Which cell classes expand or shrink with developmental age, consistently in two "
        "independent donor sets (v2 and v3 chemistry)?",
        ["Replicate unit: donor. Cells summed per donor within a tissue scope; donors with "
         f"< {C.MIN_DONOR_CELLS} cells in the scope dropped; scopes need >= {MIN_DONORS} donors "
         "in each chemistry.",
         f"Classes below {MIN_CLASS_FRACTION:.1%} of the scope merged into 'Other'; centred "
         "log-ratio with pseudocount 0.5.",
         "Spearman rho of CLR vs age per chemistry; two-sided permutation p, exact over all "
         "donor orderings (ties in age handled by permuting the observed ages).",
         "v2 and v3 combined by signed Stouffer (weights sqrt(n donors)); BH within each "
         "dataset x scope. Replicated = same direction, combined q < 0.05 and each "
         "chemistry nominally significant (one-sided p < 0.05); supported = combined "
         "q < 0.05 but one chemistry weak.",
         "Sex check: Spearman of CLR vs male indicator, reported next to each trend."],
        f,
        ["v2 spans ~6-10 pcw and v3 ~5-14 pcw with a gap at 7-11.5, so 'replicated' means "
         "monotonic across two different windows -- a rise-then-fall would fail.",
         "Composition depends on how each sample was dissected and dissociated; a trend can "
         "reflect changing dissection practice with age, which donor-level data cannot rule out.",
         "cortex and human_dev share donors, so their agreement is not independent evidence.",
         "Spearman tests monotonic trends only; n = 4-9 donors per chemistry limits power."],
        ["Add donors at the ages one chemistry lacks, so trends are tested on one window.",
         "Use per-sample dissection metadata to model composition with dissection as a covariate."])


def figures(out: C.Output, fracs: pd.DataFrame) -> None:
    plt = C.plt_or_none()
    if plt is None or fracs.empty:
        return
    for (ds, scope), g in fracs.groupby(["dataset", "scope"], sort=False):
        classes = sorted(g.cell_class.unique())
        ncol = min(4, len(classes))
        nrow = int(np.ceil(len(classes) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.6 * nrow), squeeze=False)
        for ax, cls in zip(axes.flat, classes):
            for c, marker in zip(C.CHEMISTRIES, ["o", "^"]):
                h = g[(g.cell_class == cls) & (g.chemistry == c)]
                col = np.where(h.sex == "M", "#2b6cb0", "#c05621")
                ax.scatter(h.age_pcw, h.fraction, c=col, marker=marker, s=28,
                           edgecolor="none", label=c)
            ax.set_title(cls, fontsize=9)
            ax.set_xlabel("age (pcw)", fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in list(axes.flat)[len(classes):]:
            ax.axis("off")
        axes.flat[0].set_ylabel("fraction of donor's cells", fontsize=8)
        fig.suptitle(f"{ds} / {scope}: circles v2, triangles v3; blue male, orange female",
                     fontsize=9)
        fig.tight_layout()
        name = f"composition_{ds}_{scope.replace(' ', '_')}"
        out.figure(fig, name, f"Cell-class fraction per donor vs age, {ds} / {scope}")
        plt.close(fig)


if __name__ == "__main__":
    main()
