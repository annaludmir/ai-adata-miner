#!/usr/bin/env python3
"""25 - Fit Spectra: gene programmes guided by prior knowledge, per stratum.

Spectra (Kunes et al., Nat Biotechnol 2023; package scSpectra) factorises
single-cell expression into gene programmes, steering each towards a prior gene
set through a gene-gene graph built from the sets. It keeps programmes the data
support, adapts their membership, and adds data-driven ones. Global sets apply
to all cells; cell-type sets to one class.

Prior (built here):
  global     user gene lists (config.GENE_LISTS_DIR), seed NDD and cell-cycle
             panels, and core developmental processes from GO BP (CORE_GO:
             synapse, axon guidance, Notch / Wnt / Hedgehog / BMP / FGF / RA /
             Hippo signalling, cilium, ECM, hypoxia, unfolded protein,
             glycolysis, cholesterol synthesis, translation, oxidative
             phosphorylation, splicing, chromatin remodelling, migration,
             myelination, inflammation, interferon, apoptosis); plus
             N_NEW_GLOBAL extra factors with no prior, free to find new programmes
  per class  the seed marker panel of that class (radial glia, IPC, ...)
Genes are matched to the subsample's symbols (case-insensitive, or by Ensembl
id); sets need >= MIN_SET_GENES genes detected in the subsample.

Reads <AIM_WORK>/spectra/<ns>/input.h5ad (script 24). Backend: 'cpu' (the
published full-batch model, default) or 'gpu' (the package's minibatched
Spectra_gpu module, marked by its authors as still in development).

Outputs
  <AIM_WORK>/spectra/<ns>/{model.pt, cell_scores.npz}
  csv_exports/<ns>/25_spectra/
    prior_sets.csv              every prior set: scope, genes found, genes in the model
    factor_info.csv             per factor: scope, best-matching prior set and overlap coefficient
                                of its top genes, label ('new' below OVERLAP_MIN), top genes
    factor_gene_weights.csv     factors x model genes (wide)
    factor_scores_by_group.csv  mean cell score per factor x (class | class x age | donor)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("TQDM_MININTERVAL", "60")     # Spectra's progress bar: one line a minute in job logs

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import Manifest, log
from lib.panels import load_panels

SCRIPT = "25_spectra_fit"
SUBDIR = "25_spectra"
CORE_GO = ["GO:0007411", "GO:0007416", "GO:0007268", "GO:0007219", "GO:0060070", "GO:0007224", "GO:0030509",
           "GO:0008543", "GO:0048384", "GO:0035329", "GO:0060271", "GO:0030198", "GO:0071456", "GO:0006986",
           "GO:0006096", "GO:0006695", "GO:0002181", "GO:0006120", "GO:0000398", "GO:0006338", "GO:0001764",
           "GO:0048813", "GO:0045664", "GO:0048709", "GO:0042552", "GO:0006954", "GO:0034340", "GO:0006915",
           "GO:0002040"]
MARKER_PANEL_TO_CLASS = {
    "radial_glia": "Radial glia", "neuronal_ipc": "Neuronal IPC", "neuroblast": "Neuroblast", "neuron": "Neuron",
    "glioblast_opc": "Glioblast", "oligo": "Oligo", "immune_microglia": "Immune",
    "vascular_endothelial": "Vascular", "erythrocyte": "Erythrocyte", "neural_crest": "Neural crest"}
MIN_SET_GENES = 5
N_NEW_GLOBAL = 5
OVERLAP_MIN = 0.2
N_TOP = 50


def prior(a) -> tuple[dict, pd.DataFrame]:
    upper = {g.upper(): g for g in a.var_names}
    acc = {str(x).split(".")[0]: g for x, g in zip(a.var["accession"], a.var_names)} if "accession" in a.var else {}

    def mapped(genes):
        out = []
        for g in genes:
            g = str(g).strip()
            s = upper.get(g.upper()) or acc.get(g.split(".")[0])
            if s:
                out.append(s)
        return sorted(set(out))

    rows, glob = [], {}
    panels = load_panels()
    for grp, prefix in (("user_lists", "list"), ("ndd", "seed"), ("cell_cycle", "seed")):
        for name, genes in panels.get(grp, {}).items():
            m = mapped(genes)
            rows.append({"scope": "global", "set": f"{prefix}:{name}", "n_input": len(genes), "n_found": len(m)})
            if len(m) >= MIN_SET_GENES:
                glob[f"{prefix}:{name}"] = m
    gmt = Path(config.ANNOTATIONS_DIR) / "gmt" / "GO_Biological_Process_2023.gmt"
    if gmt.exists():
        for line in gmt.read_text().splitlines():
            parts = line.split("\t")
            go = next((g for g in CORE_GO if f"({g})" in parts[0]), None)
            if go:
                m = mapped(parts[2:])
                name = "go:" + parts[0].split(" (GO:")[0]
                rows.append({"scope": "global", "set": name, "n_input": len(parts) - 2, "n_found": len(m)})
                if len(m) >= MIN_SET_GENES:
                    glob[name] = m
    else:
        log(f"  {gmt} missing -- core GO processes left out of the prior (run fetch_annotations.sh)")
    gsd = {"global": glob}
    classes = sorted(set(a.obs["cell_class"].astype(str)))
    markers = panels.get("marker", {})
    for c in classes:
        gsd[c] = {}
        for panel, cls in MARKER_PANEL_TO_CLASS.items():
            if cls == c and panel in markers:
                m = mapped(markers[panel])
                rows.append({"scope": c, "set": f"marker:{panel}", "n_input": len(markers[panel]), "n_found": len(m)})
                if len(m) >= MIN_SET_GENES:
                    gsd[c][f"marker:{panel}"] = m
    return gsd, pd.DataFrame(rows)


def run(key: str, args, chem: str | None = None, ns: str | None = None) -> None:
    import anndata as ad
    import torch
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    man = Manifest(ns, SCRIPT)
    work = config.WORK_DIR / "spectra" / ns
    if not (work / "input.h5ad").exists():
        log(f"  {work / 'input.h5ad'} missing -- run script 24 first; skipping")
        man.flush()
        return
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", torch.get_num_threads())))
    torch.manual_seed(config.RANDOM_SEED)
    np.random.seed(config.RANDOM_SEED)
    a = ad.read_h5ad(work / "input.h5ad")
    a.obs["cell_class"] = a.obs["cell_class"].astype(str)
    gsd, sets = prior(a)
    # Spectra (0.2.1) indexes .X with a pandas mask, which current scipy rejects for sparse
    # matrices, and densifies it anyway: hand it a dense matrix restricted to the genes it uses
    # (highly variable + every prior-set gene).
    in_sets = {g for v in gsd.values() for genes in v.values() for g in genes}
    use = a.var["highly_variable"].to_numpy(bool) | a.var_names.isin(list(in_sets))
    a = a[:, use].copy()
    a.X = np.asarray(a.X.todense() if hasattr(a.X, "todense") else a.X, dtype=np.float32)
    a.var["highly_variable"] = a.var["highly_variable"].astype(bool)
    L = {"global": len(gsd["global"]) + N_NEW_GLOBAL}
    for c in gsd:
        if c != "global":
            L[c] = len(gsd[c]) + 1
    log(f"  prior: {len(gsd['global'])} global sets (+{N_NEW_GLOBAL} free factors), cell-type sets "
        + ", ".join(f"{c} {len(v)}" for c, v in gsd.items() if c != "global") + f"; {sum(L.values())} factors")
    log(f"  backend {args.backend}, {args.epochs} epochs, torch threads {torch.get_num_threads()}, "
        f"cuda {'available' if torch.cuda.is_available() else 'not available'}")
    if args.backend == "gpu":
        from Spectra import Spectra_gpu as S
        model = S.est_spectra(adata=a, gene_set_dictionary=gsd, L=L, use_highly_variable=True,
                              cell_type_key="cell_class", use_weights=True, lam=args.lam, use_cell_types=True,
                              n_top_vals=N_TOP, filter_sets=True, num_epochs=args.epochs, batch_size=args.batch_size)
    else:
        from Spectra import Spectra as S
        model = S.est_spectra(adata=a, gene_set_dictionary=gsd, L=L, use_highly_variable=True,
                              cell_type_key="cell_class", use_weights=True, lam=args.lam, use_cell_types=True,
                              n_top_vals=N_TOP, filter_sets=True, label_factors=False, num_epochs=args.epochs)
    factors = np.asarray(a.uns["SPECTRA_factors"])
    scores = np.asarray(a.obsm["SPECTRA_cell_scores"])
    vocab = np.asarray(a.var_names[a.var["spectra_vocab"].to_numpy(bool)])
    im = model.internal_model
    order = next((list(v) for v in (getattr(im, "ct_order", None), getattr(im, "cell_types", None))
                  if v is not None), sorted(set(a.obs.cell_class)))
    scope = ["global"] * im.L["global"] + [c for c in order if c != "global" for _ in range(im.L[c])]
    if len(scope) != factors.shape[0]:
        log(f"  WARNING: factor scopes ({len(scope)}) do not match factors ({factors.shape[0]}); scopes left blank")
        scope = [""] * factors.shape[0]
    try:
        model.save(str(work / "model.pt"))
    except Exception as err:          # saving is a convenience; results are in the CSVs
        log(f"  model not saved: {err}")
    np.savez_compressed(work / "cell_scores.npz", scores=scores, cells=a.obs_names.to_numpy())

    flat = {s: set(g) for v in gsd.values() for s, g in v.items()}
    info = []
    for k in range(factors.shape[0]):
        top = vocab[np.argsort(-factors[k])[:N_TOP]]
        best, coef = "", 0.0
        for s, g in flat.items():
            c = len(set(top) & g) / min(len(g), N_TOP)
            if c > coef:
                best, coef = s, c
        info.append({"factor": f"F{k:02d}", "scope": scope[k], "best_prior_set": best, "overlap_coefficient": coef,
                     "label": best if coef >= OVERLAP_MIN else "new", "top_genes": "|".join(top[:20]),
                     "mean_cell_score": float(scores[:, k].mean())})
    info = pd.DataFrame(info)
    sets["in_prior"] = sets.set.isin(flat)
    man.write(sets, "prior_sets", "Prior gene sets: scope, genes in the input, genes found in the subsample",
              subdir=SUBDIR)
    man.write(info, "factor_info", f"Per factor: scope, best prior set by overlap coefficient of its top {N_TOP} "
              f"genes, label ('new' below {OVERLAP_MIN}), top genes", subdir=SUBDIR)
    w = pd.DataFrame(factors, index=info.factor, columns=vocab)
    man.write(w.reset_index(), "factor_gene_weights", "Gene weights per factor (factors x model genes)",
              subdir=SUBDIR)
    obs = a.obs.copy()
    S = pd.DataFrame(scores, index=obs.index, columns=info.factor)
    groups = {"cell_class": obs.cell_class.astype(str)}
    if "age_pcw" in obs:
        groups["cell_class_x_age"] = obs.cell_class.astype(str) + " | " + obs.age_pcw.astype(str)
    if "donor" in obs:
        groups["donor"] = obs.donor.astype(str)
    long = []
    for gname, lab in groups.items():
        m = S.groupby(lab.to_numpy()).mean()
        n = lab.value_counts()
        t = m.reset_index(names="group").melt(id_vars="group", var_name="factor", value_name="mean_score")
        t.insert(0, "grouping", gname)
        t["n_cells"] = t.group.map(n).astype(int)
        long.append(t)
    man.write(pd.concat(long, ignore_index=True), "factor_scores_by_group",
              "Mean cell score per factor per class, class x age and donor", subdir=SUBDIR)
    man.flush()


def main() -> None:
    p = cli.build_parser(__doc__)
    p.add_argument("--backend", choices=["cpu", "gpu"], default=os.environ.get("SPECTRA_BACKEND", "cpu"))
    p.add_argument("--epochs", type=int, default=int(os.environ.get("SPECTRA_EPOCHS", 0)) or None,
                   help="training epochs (default 5000 for cpu, 50 minibatch epochs for gpu)")
    p.add_argument("--batch-size", type=int, default=1000, help="gpu backend: cells per minibatch")
    p.add_argument("--lam", type=float, default=0.01, help="weight of the prior graph vs expression")
    args = p.parse_args()
    if args.epochs is None:
        args.epochs = 5000 if args.backend == "cpu" else 50
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
