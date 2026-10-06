#!/usr/bin/env python3
"""27 - Cell states within classes over development: tRG, neuron sub-types, glial precursors.

Question: 09 found that the G2/M share of cycling ventricular radial glia
(vRG) rises with age, and the oRG / vRG split did not explain it. Do truncated
radial glia (tRG, CRYAB-high, derived from vRG) appear among vRG with age, and
does the G2/M rise hold without them? Among neurons, does the balance of
deep- vs upper-layer and excitatory vs inhibitory identity shift with age; and
among glioblasts, OPC vs astrocyte-precursor identity?

Method
  A. tRG (stage-2 script 19): per stratum, the share of tRG-like cells among
     radial glia and among vRG per age (>= 50 radial glia); vRG G2/M share of
     cycling cells with and without tRG-like cells; tRG-like vs other vRG G2/M
     share at the same ages (sign test). Spearman with age, exact permutation,
     v2 x v3 combined, tiered.
  B. States (stage-2 script 22): per-cell contrasts deep - upper layer and
     excitatory - inhibitory (neurons, neuroblasts) and OPC - astrocyte
     precursor (glioblasts), each a difference of two module scores. Per class x
     age (>= 50 cells), the share of cells beyond +MARGIN (first state) and
     beyond -MARGIN (second state), and the mean contrast; Spearman with age,
     exact permutation, v2 x v3 combined, tiered.

Inputs (csv_exports/<ds>__<chem>/):
  19_rg_subtypes/{rg_trg_by_age,rg_phase_by_state_x_age}.csv
  22_cell_programs/{program_score_summary,program_score_hist}.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import _common as C

SLUG = "27_cell_states"
TITLE = "Cell states over development: truncated radial glia, neuron sub-types, glial precursors"
MIN_CELLS = 50
MIN_AGES = 5
MARGIN = 0.5
CONTRASTS = {"contrast:deep_vs_upper": (["Neuron", "Neuroblast"], "deep layer", "upper layer"),
             "contrast:excitatory_vs_inhibitory": (["Neuron", "Neuroblast"], "excitatory", "inhibitory"),
             "contrast:opc_vs_astrocyte": (["Glioblast"], "OPC", "astrocyte precursor")}


def trend(values: pd.Series) -> tuple[float, float, int]:
    v = values.dropna()
    if len(v) < MIN_AGES:
        return np.nan, np.nan, len(v)
    ages = v.index.to_numpy(float)
    rho = float(C.spearman_rows(v.to_numpy()[None, :], ages)[0])
    p, _ = C.spearman_perm_p(np.array([rho]), ages)
    return rho, float(p[0]), len(v)


def part_trg(out: C.Output) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, cmp_rows = [], []
    for ds, chem in C.STRATA:
        base = C.EXPORTS / C.ns(ds, chem) / "19_rg_subtypes"
        if not (base / "rg_trg_by_age.csv").exists():
            continue
        out.used(f"{C.ns(ds, chem)}/19_rg_subtypes/rg_trg_by_age.csv")
        t = pd.read_csv(base / "rg_trg_by_age.csv")
        t = t[(t.n_rg >= MIN_CELLS) & ~t.age.map(lambda a: C.excluded(ds, age=a))].set_index("age")
        metrics = {"tRG-like share of radial glia": t.frac_trg_like,
                   "tRG-like share of vRG": t.frac_trg_like_among_vrg.where(t.n_vrg >= 20)}
        if (base / "rg_phase_by_state_x_age.csv").exists():
            out.used(f"{C.ns(ds, chem)}/19_rg_subtypes/rg_phase_by_state_x_age.csv")
            ph = pd.read_csv(base / "rg_phase_by_state_x_age.csv")
            ph = ph[~ph.age.map(lambda a: C.excluded(ds, age=a))]
            for c in ("G1", "S", "G2M"):
                if c not in ph:
                    ph[c] = 0.0
            ph["cyc_cells"] = ph[["G1", "S", "G2M"]].sum(axis=1) * ph.n_cells
            ph["g2m_cells"] = ph["G2M"] * ph.n_cells

            def g2m_share(states):
                g = ph[ph.rg_state.isin(states)].groupby("age")[["cyc_cells", "g2m_cells"]].sum()
                return (g.g2m_cells / g.cyc_cells).where(g.cyc_cells >= 20)
            metrics["vRG G2M share of cycling (all vRG)"] = g2m_share(["vRG", "vRG tRG-like"])
            metrics["vRG G2M share of cycling (tRG-like removed)"] = g2m_share(["vRG"])
            a, b = g2m_share(["vRG tRG-like"]), g2m_share(["vRG"])
            both = pd.concat([a, b], axis=1, keys=["trg", "other"]).dropna()
            if len(both):
                k = int((both.trg > both.other).sum())
                cmp_rows.append({"dataset": ds, "chemistry": chem, "n_ages": len(both),
                                 "trg_higher_g2m_at": k, "mean_trg": float(both.trg.mean()),
                                 "mean_other_vrg": float(both.other.mean()),
                                 "sign_p": float(stats.binomtest(k, len(both)).pvalue)})
        for m, v in metrics.items():
            rho, p, k = trend(v)
            if np.isfinite(rho):
                vv = v.dropna()
                rows.append({"dataset": ds, "chemistry": chem, "metric": m, "n_ages": k, "rho_vs_age": rho,
                             "perm_p": p, "range": f"{vv.iloc[0]:.3g}->{vv.iloc[-1]:.3g}"})
    return pd.DataFrame(rows), pd.DataFrame(cmp_rows)


def part_states(out: C.Output) -> pd.DataFrame:
    rows = []
    for ds, chem in C.STRATA:
        base = C.EXPORTS / C.ns(ds, chem) / "22_cell_programs"
        if not (base / "program_score_hist.csv").exists():
            continue
        out.used(f"{C.ns(ds, chem)}/22_cell_programs/program_score_hist.csv",
                 f"{C.ns(ds, chem)}/22_cell_programs/program_score_summary.csv")
        h = pd.read_csv(base / "program_score_hist.csv")
        s = pd.read_csv(base / "program_score_summary.csv")
        h = h[(h.grouping == "cell_class_x_age") & h.program.isin(CONTRASTS)]
        s = s[(s.grouping == "cell_class_x_age") & s.program.isin(CONTRASTS)].set_index(["group", "program"])
        width = 0.05
        for (grp, prog), g in h.groupby(["group", "program"]):
            cls, age = C.parse_group(grp, 2)
            classes, a_name, b_name = CONTRASTS[prog]
            try:
                age = float(age)
            except ValueError:
                continue
            if cls not in classes or C.excluded(ds, age=age):
                continue
            n = g["count"].sum()
            if n < MIN_CELLS:
                continue
            lo = g.bin_lo.to_numpy(float)
            first = g.loc[lo >= MARGIN, "count"].sum() / n
            second = g.loc[(lo + width) <= -MARGIN, "count"].sum() / n
            rows.append({"dataset": ds, "chemistry": chem, "contrast": prog.split(":")[1], "cell_class": cls,
                         "age": age, "n_cells": int(n), "share_first_state": first, "share_second_state": second,
                         "first_state": a_name, "second_state": b_name,
                         "mean_contrast": float(s.loc[(grp, prog), "mean"]) if (grp, prog) in s.index else np.nan})
    return pd.DataFrame(rows)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    trg, trg_cmp = part_trg(out)
    states = part_states(out)
    if trg.empty and states.empty:
        out.summary(TITLE, "Skipped: no tRG flags or state scores yet.",
                    ["Needs 19_rg_subtypes/rg_trg_by_age.csv and 22_cell_programs/ from stage 2."],
                    ["**Not run**: the tRG flag (script 19) and state scores (script 22) are new; re-run stage 2 "
                     "(slurm_02_pseudobulk.sh), then this analysis."], [], [])
        return
    tcomb = C.combine_chemistries(trg, ["dataset", "metric"], effect="rho_vs_age", weight="n_ages",
                                  labels=("rises with age", "falls with age"), carry=("rho_vs_age", "range")) \
        if len(trg) else pd.DataFrame()
    srows = []
    for (ds, chem, con, cls), g in states.groupby(["dataset", "chemistry", "contrast", "cell_class"]) if len(states) else []:
        g = g.set_index("age").sort_index()
        for col, label in (("share_first_state", f"share {g.first_state.iloc[0]}"),
                           ("share_second_state", f"share {g.second_state.iloc[0]}"),
                           ("mean_contrast", "mean contrast")):
            rho, p, k = trend(g[col])
            if np.isfinite(rho):
                vv = g[col].dropna()
                srows.append({"dataset": ds, "chemistry": chem, "contrast": con, "cell_class": cls, "metric": label,
                              "n_ages": k, "rho_vs_age": rho, "perm_p": p, "range": f"{vv.iloc[0]:.3g}->{vv.iloc[-1]:.3g}"})
    sper = pd.DataFrame(srows)
    scomb = C.combine_chemistries(sper, ["dataset", "contrast", "cell_class", "metric"], effect="rho_vs_age",
                                  weight="n_ages", labels=("rises with age", "falls with age"),
                                  carry=("rho_vs_age", "range")) if len(sper) else pd.DataFrame()
    out.write(trg, "trg_trends_per_stratum", "Per stratum: tRG-like shares and vRG G2/M share vs age")
    out.write(tcomb, "trg_trends_combined", "v2 x v3 combined; tier")
    out.write(trg_cmp, "trg_vs_other_vrg_g2m", "G2/M share of cycling tRG-like vs other vRG at the same ages (sign test)")
    out.write(states, "state_shares_by_age", f"Per stratum x class x age: share of cells beyond +/-{MARGIN} on each contrast")
    out.write(sper, "state_trends_per_stratum", "Per stratum x contrast x class: state shares vs age")
    out.write(scomb, "state_trends_combined", "v2 x v3 combined; tier")

    f = []
    if len(tcomb):
        f.append("**tRG and the vRG G2/M rise** (rho v2/v3; range v2 / v3): " + "; ".join(
            f"{r.dataset} {r.metric} {r.direction} ({r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f}; "
            f"{r.range_v2} / {r.range_v3}; {r.tier or 'n.s.'})" for r in tcomb.itertuples()) + ".")
        w = tcomb.set_index(["dataset", "metric"])
        for ds in C.DATASETS:
            a, b = (ds, "vRG G2M share of cycling (all vRG)"), (ds, "vRG G2M share of cycling (tRG-like removed)")
            if a in w.index and b in w.index:
                f.append(f"**{ds}: does the vRG G2/M trend survive removing tRG-like cells?** all vRG "
                         f"{w.loc[a, 'tier'] or 'n.s.'} -> tRG-like removed {w.loc[b, 'tier'] or 'n.s.'}"
                         + (" -- tRG do not explain it." if w.loc[a, "tier"] and w.loc[b, "tier"] else "."))
    if len(trg_cmp):
        f.append("**tRG-like vs other vRG, G2/M share at the same ages**: " + "; ".join(
            f"{r.dataset} {r.chemistry} {r.mean_trg:.2f} vs {r.mean_other_vrg:.2f}, tRG-like higher at "
            f"{r.trg_higher_g2m_at}/{r.n_ages} ages (sign p = {C.fmt_p(r.sign_p)})" for r in trg_cmp.itertuples()) + ".")
    if len(scomb):
        rep = scomb[scomb.tier != ""]
        f.append("**State balance changing with age** (rho v2/v3; range v2 / v3): " + (
            "; ".join(f"{r.dataset} {r.cell_class} {r.metric} ({r.contrast.replace('_', ' ')}) {r.direction} "
                      f"({r.rho_vs_age_v2:+.2f}/{r.rho_vs_age_v3:+.2f}; {r.range_v2} / {r.range_v3}; {r.tier})"
                      for r in rep.itertuples()) if len(rep) else "none replicated") + ".")
    out.summary(
        TITLE,
        "Do tRG appear among vRG with age, and does the vRG G2/M rise hold without them? Does the balance of "
        "neuron sub-types and of glial-precursor states shift with age?",
        ["A: tRG-like = non-oRG radial glia with CRYAB >= 1.5 log1p(CP10K) (script 19); shares per age (>= "
         f"{MIN_CELLS} radial glia); vRG G2/M share with and without tRG-like cells; Spearman with age, exact "
         "permutation, v2 x v3 signed Stouffer, tiered.",
         f"B: per-cell contrasts of module scores (script 22); share of cells beyond +/-{MARGIN} log units per "
         f"class x age (>= {MIN_CELLS} cells); trends as in A."],
        f,
        ["tRG are reported from mid-gestation (around 16-17 pcw), later than most ages here; a CRYAB flag at "
         "earlier ages marks CRYAB-expressing radial glia, not necessarily tRG.",
         "Contrasts of short marker programmes are noisy per cell; the shares are best read as trends.",
         "Deep-layer neurons are born first, so a deep-to-upper shift with age is expected and serves as a check."],
        ["Score cells with a reference atlas (e.g. label transfer) instead of marker programmes."])


if __name__ == "__main__":
    main()
