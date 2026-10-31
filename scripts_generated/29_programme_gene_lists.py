#!/usr/bin/env python3
"""29 - Export gene programmes as plain gene lists (for GSEA in ndd_gene_modules).

Reads the programme fits' CSV exports (no h5ad) and writes one CSV per
programme with its N_TOP highest-weighted genes, in the format the
ndd_gene_modules gene-list tools read:

    gene
    TOP2A
    CDK1
    ...

Default fit: the core-prior Spectra fit (25_spectra_core), whose prior holds
GO processes, cell-cycle panels and class markers but no NDD list or NDD
panel. Its programmes are of two kinds, kept apart in the file names and the
index: 'free' ones (no prior set at all, learned from the data alone) and
'steered' ones (guided by a core-biology set). --free-only keeps the free ones.
--fit string / cnmf / schpf export the STRING-prior Spectra fit or the
prior-free fits instead.

Replication: v2 and v3 programmes of a dataset are paired when each is the
other's best match by cosine of gene weights (>= MATCH_MIN), as in analysis 36.
The index records each programme's partner; --replicated-only keeps paired
programmes only, and for each pair a consensus list is also written (genes in
the top N of both, ordered by mean rank).

Outputs (csv_exports/_programme_gene_lists/<fit>/; replaced on every run)
  <dataset>__<chem>__<factor>__<kind>.csv   one gene list per programme
  consensus/<dataset>__<f_v2>~<f_v3>.csv    genes shared by a replicated pair
  index.csv                                 per list: dataset, chemistry, factor, kind, scope,
                                            label, partner, cosine, n_genes, top genes, file
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import log

SCRIPT = "29_programme_gene_lists"
FITS = {"core": "25_spectra_core", "string": "25_spectra_string", "full": "25_spectra",
        "cnmf": "26_cnmf", "schpf": "27_schpf"}
N_TOP = 100
MATCH_MIN = 0.5


def load_fit(ns: str, sub: str):
    d = config.CSV_EXPORTS / ns / sub
    if not ((d / "factor_info.csv").exists() and (d / "factor_gene_weights.csv").exists()):
        return None
    info = pd.read_csv(d / "factor_info.csv").set_index("factor")
    w = pd.read_csv(d / "factor_gene_weights.csv").set_index("factor")
    if "free" not in info:              # older fits: free = last of each scope (5 global / 1 per class)
        info["free"] = False
        for scope, g in info.groupby("scope", sort=False):
            info.loc[g.index[-(5 if scope == "global" else 1):], "free"] = True
    return info, w


def pair(w2: pd.DataFrame, w3: pd.DataFrame) -> dict:
    genes = w2.columns.intersection(w3.columns)
    a, b = w2[genes].to_numpy(float), w3[genes].to_numpy(float)
    a = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
    b = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-12)
    S = a @ b.T
    best3, best2 = S.argmax(axis=1), S.argmax(axis=0)
    out = {}
    for i, j in enumerate(best3):
        if best2[j] == i and S[i, j] >= MATCH_MIN:
            out[("v2", w2.index[i])] = ("v3", w3.index[j], float(S[i, j]))
            out[("v3", w3.index[j])] = ("v2", w2.index[i], float(S[i, j]))
    return out


def safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")


def main() -> None:
    p = cli.build_parser(__doc__)
    p.add_argument("--fit", choices=list(FITS), default="core", help="which fit to export (default: core prior)")
    p.add_argument("--n-top", type=int, default=N_TOP, help="genes per programme")
    p.add_argument("--free-only", action="store_true", help="only programmes learned without any prior set")
    p.add_argument("--replicated-only", action="store_true", help="only programmes paired between v2 and v3")
    args = p.parse_args()
    sub = FITS[args.fit]
    out = config.CSV_EXPORTS / "_programme_gene_lists" / args.fit
    (out / "consensus").mkdir(parents=True, exist_ok=True)
    for old in list(out.glob("*.csv")) + list((out / "consensus").glob("*.csv")):
        old.unlink()                      # lists of an earlier run with other options

    rows, cons = [], []
    for ds in cli.selected_datasets(args):
        fits = {ch: load_fit(config.namespace(ds, ch), sub) for ch in ("v2", "v3")}
        if all(f is None for f in fits.values()):
            log(f"{ds}: no {sub} fits -- skipping")
            continue
        pairs = pair(fits["v2"][1], fits["v3"][1]) if all(f is not None for f in fits.values()) else {}
        tops = {}
        for ch, fit in fits.items():
            if fit is None:
                continue
            info, w = fit
            for f in w.index:
                free = bool(info.loc[f, "free"])
                partner = pairs.get((ch, f))
                if (args.free_only and not free) or (args.replicated_only and partner is None):
                    continue
                genes = list(w.loc[f].sort_values(ascending=False).index[:args.n_top])
                tops[(ch, f)] = genes
                kind = "free" if free else "steered-" + safe(info.loc[f, "label"])
                name = f"{ds}__{ch}__{f}__{kind}.csv"
                pd.DataFrame({"gene": genes}).to_csv(out / name, index=False)
                rows.append({"dataset": ds, "chemistry": ch, "factor": f, "kind": "free" if free else "steered",
                             "scope": info.loc[f, "scope"], "label": info.loc[f, "label"],
                             "partner": f"{partner[0]} {partner[1]}" if partner else "",
                             "cosine": partner[2] if partner else np.nan, "n_genes": len(genes),
                             "top_genes": "|".join(genes[:10]), "file": name})
        for (ch, f), (_, g, cos) in pairs.items():
            if ch != "v2" or ("v2", f) not in tops or ("v3", g) not in tops:
                continue
            r2 = {x: i for i, x in enumerate(tops[("v2", f)])}
            r3 = {x: i for i, x in enumerate(tops[("v3", g)])}
            shared = sorted(set(r2) & set(r3), key=lambda x: r2[x] + r3[x])
            name = f"{ds}__{f}~{g}.csv"
            pd.DataFrame({"gene": shared}).to_csv(out / "consensus" / name, index=False)
            cons.append({"dataset": ds, "pair": f"{f}~{g}", "cosine": cos, "n_shared": len(shared),
                         "file": f"consensus/{name}"})
        log(f"{ds}: {sum(r['dataset'] == ds for r in rows)} programme lists, "
            f"{sum(c['dataset'] == ds for c in cons)} replicated-pair consensus lists")
    pd.DataFrame(rows).to_csv(out / "index.csv", index=False)
    pd.DataFrame(cons).to_csv(out / "consensus" / "index.csv", index=False)
    log(f"wrote {len(rows)} lists to {out}")


if __name__ == "__main__":
    main()
