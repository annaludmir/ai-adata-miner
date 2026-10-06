#!/usr/bin/env python3
"""20 - Do the files' other cell-type annotations agree with CellClass?

Question: every analysis here groups cells by CellClass. Each file carries a
second annotation -- cortex an older label ('classes'), human_dev Cell
Ontology terms (cell_type_ontology_term_id). Where do they agree, which
CellClass groups does each other label split or merge, and is the mapping the
same in both chemistries? Disagreement marks the classes whose boundaries
deserve caution (for example IPC vs neuroblast).

Method
  Cross-tabulations CellClass x other annotation from stage 1 (02_composition,
  roles alt_cell_class and cell_type_id; exported from this version of the
  repo on). Per stratum: adjusted Rand index and normalised mutual
  information between the two labelings; for each other-label level, its
  dominant CellClass and that class's share; for each CellClass, its dominant
  other label. Mappings compared between v2 and v3.

Inputs (csv_exports/):
  <ds>__<chem>/02_composition/counts_cell_class_by_{alt_cell_class,cell_type_id}.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _common as C

SLUG = "20_annotation_agreement"
TITLE = "Agreement between CellClass and the files' other cell-type annotations"
ALT_ROLES = {"alt_cell_class": "older 'classes' label", "cell_type_id": "Cell Ontology term"}
# Cell Ontology ids seen in the atlas, named where the mapping is certain
CL_NAMES = {"CL:0000540": "neuron", "CL:0000681": "radial glial cell", "CL:0000031": "neuroblast",
            "CL:0000030": "glioblast", "CL:0000128": "oligodendrocyte",
            "CL:0002453": "oligodendrocyte precursor cell", "CL:0000129": "microglial cell",
            "CL:0000115": "endothelial cell", "CL:0000232": "erythrocyte", "CL:0000057": "fibroblast",
            "CL:0000669": "pericyte", "CL:0000127": "astrocyte", "CL:0000738": "leukocyte"}


def comb2(x: np.ndarray) -> np.ndarray:
    return x * (x - 1) / 2


def ari(ct: np.ndarray) -> float:
    n = ct.sum()
    a, b = comb2(ct.sum(axis=1)).sum(), comb2(ct.sum(axis=0)).sum()
    expected = a * b / comb2(n)
    top = (a + b) / 2 - expected
    return float((comb2(ct).sum() - expected) / top) if top > 0 else np.nan


def nmi(ct: np.ndarray) -> float:
    p = ct / ct.sum()
    pi, pj = p.sum(axis=1), p.sum(axis=0)
    nz = p > 0
    mi = (p[nz] * np.log(p[nz] / np.outer(pi, pj)[nz])).sum()
    h = lambda q: -(q[q > 0] * np.log(q[q > 0])).sum()
    den = np.sqrt(h(pi) * h(pj))
    return float(mi / den) if den > 0 else np.nan


def has_column(namespace: str, role: str) -> bool:
    """Does the file carry a column for this role (01_overview counts every obs column)?"""
    if str(C.REPO) not in sys.path:
        sys.path.insert(0, str(C.REPO))
    import config
    return any((C.EXPORTS / namespace / "01_overview" / f"cell_counts_{col}.csv").exists()
               for col in config.COLUMN_ROLES.get(role, []))


def name(level: str) -> str:
    return f"{level} ({CL_NAMES[level]})" if level in CL_NAMES else str(level)


def main() -> None:
    out = C.Output(SLUG)
    C.log(f"=== {SLUG}")
    summary, levels, missing = [], [], []
    for ds, chem in C.STRATA:
        n = C.ns(ds, chem)
        for role, label in ALT_ROLES.items():
            p = C.EXPORTS / n / "02_composition" / f"counts_cell_class_by_{role}.csv"
            if not p.exists():
                if has_column(n, role):
                    missing.append(f"{n} {role}")
                continue
            out.used(f"{n}/02_composition/counts_cell_class_by_{role}.csv")
            ct = pd.read_csv(p).set_index(role)
            ct.index = ct.index.astype(str)
            ct = ct.loc[ct.sum(axis=1) > 0, ct.sum(axis=0) > 0]
            M = ct.to_numpy(float)
            summary.append({"dataset": ds, "chemistry": chem, "annotation": role, "description": label,
                            "n_levels": ct.shape[0], "n_cell_classes": ct.shape[1], "n_cells": int(M.sum()),
                            "adjusted_rand_index": ari(M), "normalised_mutual_information": nmi(M),
                            "share_in_dominant_pairing": float(M.max(axis=1).sum() / M.sum())})
            for lvl, row in ct.iterrows():
                top = row.sort_values(ascending=False)
                levels.append({"dataset": ds, "chemistry": chem, "annotation": role, "level": lvl,
                               "level_name": name(lvl), "n_cells": int(row.sum()),
                               "dominant_cell_class": top.index[0], "dominant_share": top.iloc[0] / row.sum(),
                               "second_cell_class": top.index[1] if len(top) > 1 and top.iloc[1] > 0 else "",
                               "second_share": top.iloc[1] / row.sum() if len(top) > 1 else 0.0})
            for cls in ct.columns:
                col = ct[cls].sort_values(ascending=False)
                levels.append({"dataset": ds, "chemistry": chem, "annotation": f"{role} (per CellClass)",
                               "level": cls, "level_name": cls, "n_cells": int(col.sum()),
                               "dominant_cell_class": name(col.index[0]),
                               "dominant_share": col.iloc[0] / col.sum(),
                               "second_cell_class": name(col.index[1]) if len(col) > 1 and col.iloc[1] > 0 else "",
                               "second_share": col.iloc[1] / col.sum() if len(col) > 1 else 0.0})
    summary, levels = pd.DataFrame(summary), pd.DataFrame(levels)
    if summary.empty:
        out.summary(TITLE, "Skipped: the cross-tabulations are not exported yet.",
                    ["Needs 02_composition's counts_cell_class_by_{alt_cell_class,cell_type_id}.csv."],
                    ["**Not run**: these tables are written by stage 1 from this version of the repo on "
                     "(config role alt_cell_class, and cell_type_id added to 02_composition's groupings). "
                     "Re-run stage 1 (slurm_01_metadata.sh), then this analysis."],
                    [], [])
        return
    out.write(summary, "annotation_agreement", "Per stratum x annotation: ARI, NMI, share of cells in the "
              "dominant pairing")
    out.write(levels, "annotation_mapping", "Per annotation level: dominant CellClass and its share; and per "
              "CellClass: dominant other label")

    f = []
    for (ds, role), g in summary.groupby(["dataset", "annotation"], sort=False):
        f.append(f"**{ds}, {ALT_ROLES[role]} vs CellClass**: ARI " + "/".join(f"{x:.2f}" for x in g.adjusted_rand_index)
                 + ", NMI " + "/".join(f"{x:.2f}" for x in g.normalised_mutual_information)
                 + " (v2/v3); " + "/".join(f"{x:.0%}" for x in g.share_in_dominant_pairing)
                 + " of cells in their level's dominant class.")
        lv = levels[(levels.dataset == ds) & (levels.annotation == role) & (levels.n_cells >= 100)]
        split = lv[lv.dominant_share < 0.8].drop_duplicates("level")
        if len(split):
            f.append(f"**{ds}: {ALT_ROLES[role]} levels that straddle CellClass groups** (dominant share < 80%): "
                     + "; ".join(f"{r.level_name} = {r.dominant_cell_class} {r.dominant_share:.0%} + "
                                 f"{r.second_cell_class} {r.second_share:.0%}" for r in split.head(10).itertuples()) + ".")
        w = lv.pivot_table(index="level", columns="chemistry", values="dominant_cell_class", aggfunc="first").dropna()
        if set(C.CHEMISTRIES) <= set(w.columns):
            diff = w[w["v2"] != w["v3"]]
            f.append(f"**{ds}: mapping consistency** -- {len(w) - len(diff)} of {len(w)} {ALT_ROLES[role]} levels "
                     "map to the same dominant CellClass in v2 and v3"
                     + (": differ for " + ", ".join(name(x) for x in diff.index[:8]) if len(diff) else "") + ".")
    if missing:
        f.append("**Not exported yet** (re-run stage 1): " + ", ".join(missing) + ".")
    out.summary(
        TITLE,
        "Do the files' other cell-type annotations agree with CellClass, which classes do they split or "
        "merge, and is the mapping the same in both chemistries?",
        ["Cross-tabulations from stage 1 (02_composition); adjusted Rand index and normalised mutual "
         "information per stratum; dominant class per level and per CellClass; v2 vs v3 mapping."],
        f,
        ["Agreement between two labels is not correctness: both may share a boundary choice.",
         "Cell Ontology names are given only for ids whose mapping is certain; others are shown as ids."],
        ["Score cells with marker panels (11_panels) to adjudicate levels that straddle CellClass groups."])


if __name__ == "__main__":
    main()
