#!/usr/bin/env python3
"""Figures for docs/FINDINGS.md, drawn from the step-3 result tables.

    python docs/figures/make_figures.py [RESULTS_DIR]

RESULTS_DIR defaults to $AIM_DOWNSTREAM_OUT or ~/aim_downstream/results. Writes
fig1..fig6 PNGs next to this script. Colours: the dataviz reference palette
(v2 = categorical slot 1 blue, v3 = slot 2 orange; diverging blue - grey - red).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

RES = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("AIM_DOWNSTREAM_OUT",
                                                                  Path.home() / "aim_downstream" / "results"))
OUT = Path(__file__).resolve().parent

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
V2, V3 = "#2a78d6", "#eb6834"
DIVERGING = LinearSegmentedColormap.from_list("bgr", ["#104281", "#3987e5", "#cde2fb", "#f0efec",
                                                      "#f6c9c0", "#e34948", "#8f2422"])
CLASS_ORDER = ["Radial glia", "Neuronal IPC", "Neuroblast", "Neuron", "Glioblast", "Oligo", "Immune",
               "Vascular", "Fibroblast", "Erythrocyte", "Neural crest", "Placodes"]
LISTS = ["ASD_HC", "autism_strong", "ASD_SFARI_syndromic", "ASD_SFARI_2", "ASD_SFARI_3", "microcephaly_genes",
         "BP_GWAS", "AD_GWAS"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.titleweight": "bold",
    "axes.titlesize": 11, "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "legend.frameon": False,
})


def save(fig, name):
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


def fig1_design():
    t = pd.read_csv(RES / "01_data_audit" / "replicate_structure.csv")
    order = ["cortex__v2", "cortex__v3", "human_dev__v2", "human_dev__v3"]
    t = t.set_index("stratum").loc[order]
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    for i, (s, r) in enumerate(t.iterrows()):
        ages = [float(a) for a in str(r.ages).split("|")]
        col = V2 if s.endswith("v2") else V3
        ax.scatter(ages, [i] * len(ages), s=46, color=col, zorder=3, edgecolor=SURFACE, linewidth=1.5)
        ax.text(14.6, i, f"{int(r.n_donors)} donors, {int(r.n_cells):,} cells", va="center", fontsize=9, color=INK2)
    ax.set_yticks(range(len(order)), [s.replace("__", " ") for s in order])
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_xlim(4.5, 14.4)
    ax.set_xlabel("Age (post-conception weeks)")
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_title("Two datasets, each split by 10x chemistry into two donor-disjoint replicates")
    save(fig, "fig1_design.png")


def _heat(ax, M, R, title, vlim, fmt="{:+.1f}"):
    im = ax.imshow(M.to_numpy(float), cmap=DIVERGING, vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(M.shape[1]), M.columns, rotation=35, ha="right")
    ax.set_yticks(range(M.shape[0]), M.index)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iat[i, j]
            if np.isfinite(v):
                rep = R.iat[i, j]
                ax.text(j, i, fmt.format(v) + ("*" if rep else ""), ha="center", va="center", fontsize=7.5,
                        color="white" if abs(v) > 0.6 * vlim else INK, fontweight="bold" if rep else "normal")
    ax.set_title(title)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    return im


def fig2_where():
    t = pd.read_csv(RES / "06_gene_list_landscape" / "list_class_preference_combined.csv")
    t["mean_effect"] = t[["effect_v2", "effect_v3"]].mean(axis=1)
    t["rep"] = t.tier.fillna("") == "replicated"
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4), gridspec_kw={"width_ratios": [5, 10]})
    for ax, ds in zip(axes, ["cortex", "human_dev"]):
        g = t[(t.dataset == ds) & t.gene_list.isin(LISTS)]
        cls = [c for c in CLASS_ORDER if c in set(g.cell_class)]
        M = g.pivot(index="gene_list", columns="cell_class", values="mean_effect").reindex(index=LISTS, columns=cls)
        R = g.pivot(index="gene_list", columns="cell_class", values="rep").reindex(index=LISTS, columns=cls).fillna(False)
        im = _heat(ax, M, R, ds, 9)
        if ds == "human_dev":
            ax.set_yticks([])
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.01)
    cb.set_label("Enrichment vs expression-matched genes\n(null SDs, mean of v2 and v3)", color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle("Where each NDD gene list is expressed (* replicated in both chemistries)", x=0.01, ha="left",
                 fontweight="bold", fontsize=12)
    save(fig, "fig2_lists_by_cell_class.png")


def fig3_age():
    t = pd.read_csv(RES / "06_gene_list_landscape" / "list_age_coordination_combined.csv")
    t = t[t.dataset == "cortex"].copy()
    t["mean_rho"] = t[["mean_rho_v2", "mean_rho_v3"]].mean(axis=1)
    t["rep"] = t.tier.fillna("") == "replicated"
    cls = ["Radial glia", "Neuronal IPC", "Neuroblast", "Neuron"]
    M = t.pivot(index="gene_list", columns="cell_class", values="mean_rho").reindex(index=LISTS, columns=cls)
    R = t.pivot(index="gene_list", columns="cell_class", values="rep").reindex(index=LISTS, columns=cls).fillna(False)
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    im = _heat(ax, M, R, "", 0.5, fmt="{:+.2f}")
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("Mean gene rho with age\n(mean of v2 and v3)", color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle("Cortex: ASD lists rise with age inside progenitor and neuroblast classes", x=0.01, ha="left",
                 fontweight="bold", fontsize=12, y=1.04)
    fig.text(0.01, 0.975, "* replicated in both chemistries against expression-matched random genes",
             fontsize=8.5, color=INK2, ha="left")
    save(fig, "fig3_list_age_trends_cortex.png")


def fig4_proliferation():
    t = pd.read_csv(RES / "09_cell_cycle_programs" / "proliferation_trajectories_per_stratum.csv")
    t = t[(t.dataset == "cortex") & (t.cell_class == "Radial glia")]
    panels = [("oRG share of radial glia", "Outer radial glia share\nof radial glia"),
              ("G2M share of cycling cells (fixed sub-type mix)", "G2/M share of cycling radial glia\n(sub-type mix held fixed)"),
              ("fraction cycling", "Fraction of radial glia\nthat are cycling")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=False)
    for ax, (metric, label) in zip(axes, panels):
        for chem, col in (("v2", V2), ("v3", V3)):
            r = t[(t.chemistry == chem) & (t.metric == metric)].iloc[0]
            a0, a1 = [float(x) for x in r.age_range.split("-")]
            ax.plot([a0, a1], [r["first"], r["last"]], color=col, linewidth=2, marker="o", markersize=7,
                    markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
            ax.text(a1 + 0.25, r["last"], f"{chem}  rho {r.rho:+.2f}", color=INK2, fontsize=8.5, va="center")
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Age (pcw), youngest and oldest donor")
        ax.set_xlim(4.5, 17)
        ax.yaxis.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
    fig.suptitle("Cortex radial glia: more outer radial glia, a larger G2/M share and fewer cycling cells with age",
                 x=0.01, ha="left", fontweight="bold", fontsize=12, y=1.06)
    save(fig, "fig4_radial_glia_proliferation.png")


def fig5_programmes():
    rows = []
    for fit, name in (("core", "Spectra, core prior"), ("string", "Spectra, STRING prior"), ("schpf", "scHPF (no prior)"),
                      ("cnmf", "cNMF (no prior)")):
        r = pd.read_csv(RES / "36_spectra_programs" / f"{fit}_ndd_enrichment_replicated.csv")
        r = r[(r.dataset == "human_dev") & (r.gene_set == "list:microcephaly_genes")]
        r = r.assign(fmin=r[["fold_v2", "fold_v3"]].min(axis=1)).sort_values("fmin", ascending=False).iloc[0]
        rows.append((name, r.fold_v2, r.fold_v3, r.q_v2, r.q_v3, bool(r.replicated)))
    fig, ax = plt.subplots(figsize=(7.4, 3.3))
    y = np.arange(len(rows))
    for i, (name, f2, f3, q2, q3, rep) in enumerate(rows):
        ax.barh(i - 0.18, f2, height=0.34, color=V2)
        ax.barh(i + 0.18, f3, height=0.34, color=V3)
        ax.text(max(f2, f3) + 0.3, i, "replicated" if rep else f"q {q2:.2f} / {q3:.4f}", va="center",
                fontsize=8.5, color=INK2)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.invert_yaxis()
    ax.axvline(1, color=MUTED, linewidth=1, linestyle=(0, (3, 3)))
    ax.set_xlabel("Microcephaly genes among the programme's top 50 genes (fold over expression-matched expectation)")
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(0, 16.5)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=V2, label="v2 donors"), Patch(color=V3, label="v3 donors")], loc="lower right")
    ax.set_title("human_dev: the progenitor proliferation programme carries the microcephaly genes in every model")
    save(fig, "fig5_microcephaly_programme.png")


def fig6_brainspan():
    t = pd.read_csv(RES / "34_brainspan_validation" / "set_trends_ours_vs_brainspan.csv")
    t = t[t.gene_set.str.match(r"^(list|seed):") & (t.gene_set != "list:BP_GWAS_H3_annotation")]
    t = t.sort_values("mean_ours").reset_index(drop=True)
    same = t.same_direction.astype(str).str.lower().isin(["true", "yes", "1"])
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    y = np.arange(len(t))
    ax.axvline(0, color=AXIS, linewidth=1)
    for i, r in t.iterrows():
        ax.plot([r.mean_ours, r.mean_brainspan], [i, i], color=GRID, linewidth=2, zorder=1)
    ax.scatter(t.mean_ours, y, s=52, color=V2, edgecolor=SURFACE, linewidth=1.5, zorder=3, label="this study (cortex)")
    ax.scatter(t.mean_brainspan, y, s=52, marker="D", color="#1baf7a", edgecolor=SURFACE, linewidth=1.5, zorder=3,
               label="BrainSpan neocortex")
    ax.set_yticks(y, [g.split(":", 1)[1] + ("" if ok else "  (opposite)") for g, ok in zip(t.gene_set, same)])
    ax.set_xlabel("Mean gene rho with age")
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8.5)
    ax.set_title("Gene-set age trends point the same way in BrainSpan (8-16 pcw) for 12 of 15 sets")
    save(fig, "fig6_brainspan_validation.png")


if __name__ == "__main__":
    for f in (fig1_design, fig2_where, fig3_age, fig4_proliferation, fig5_programmes, fig6_brainspan):
        f()
