#!/usr/bin/env python3
"""05 - Which covariates are entangled with which?

The single most dangerous thing about a developmental atlas is that age,
chemistry and donor are usually nested: a "developmental" expression change can
simply be the v2 -> v3 chemistry switch.  This quantifies every such
entanglement up front so later analyses can be read with the right caveats.

  * categorical x categorical -> Cramer's V (bias-corrected)
  * categorical x numeric QC  -> eta-squared (variance explained)
  * plus the raw contingency tables for the pairs that matter most

Outputs (csv_exports/<dataset>/05_confounds/)
  covariate_association_matrix.csv  square Cramer's V matrix
  covariate_association_pairs.csv   long form, sorted, with a confound verdict
  qc_variance_explained.csv         eta-squared of each QC metric per covariate
  crosstab_<a>_x_<b>.csv            contingency tables for key pairs
  confound_warnings.csv             pairs that are strongly entangled
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.aggregate import contingency
from lib.io_utils import (Manifest, load_obs, log,
                          resolve_cluster_columns, resolve_qc_frame, resolve_role)
from lib.stats_utils import cramers_v, eta_squared

SCRIPT = "05_covariate_confounds"
SUBDIR = "05_confounds"
COVARIATE_ROLES = ["age", "chemistry", "donor", "sample", "region", "subregion",
                   "dissection", "cell_class", "cyclephase", "sex", "dev_stage",
                   "assay", "transcriptome", "batch"]
KEY_PAIRS = [("age", "chemistry"), ("age", "donor"), ("donor", "chemistry"),
             ("region", "donor"), ("region", "chemistry"), ("age", "region"),
             ("cell_class", "chemistry"), ("cell_class", "donor")]
MAX_LEVELS = 400
STRONG_V = 0.5


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)

    obs, _keep = load_obs(path, key, chem, args.limit_cells)
    if obs is None:
        log(f"  no chemistry column in this file -- cannot run chemistry={chem}; skipping")
        man.flush()
        return
    if len(obs) == 0:
        log(f"  no cells with chemistry={chem}; skipping")
        man.flush()
        return
    log(f"  {len(obs):,} cells after chemistry filter ({chem or 'pooled'})")
    qc = resolve_qc_frame(obs)

    covars: dict[str, pd.Series] = {}
    for role in COVARIATE_ROLES:
        col = "age_pcw" if role == "age" and "age_pcw" in obs.columns else resolve_role(obs, role)
        if role == "cyclephase" and "cyclephase_h" in obs.columns:
            col = "cyclephase_h"
        if col is None or col not in obs.columns:
            continue
        s = obs[col].astype("object")
        n = s.nunique(dropna=True)
        if n < 2 or n > MAX_LEVELS:
            if n > MAX_LEVELS:
                log(f"  skipping covariate '{role}' ({col}): {n:,} levels")
            continue
        covars[role] = s
    for col in resolve_cluster_columns(obs):
        if 2 <= obs[col].nunique(dropna=True) <= MAX_LEVELS:
            covars[f"cluster_{col}"] = obs[col].astype("object")

    log(f"  covariates compared: {', '.join(covars)}")
    if len(covars) < 2:
        log("  fewer than two usable covariates; nothing to compare")
        man.flush()
        return

    names = list(covars)
    pairs = []
    matrix = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for a, b in combinations(names, 2):
        tab = pd.crosstab(covars[a], covars[b])
        v = cramers_v(tab)
        matrix.loc[a, b] = matrix.loc[b, a] = v
        pairs.append({
            "covariate_a": a, "covariate_b": b,
            "n_levels_a": int(covars[a].nunique()), "n_levels_b": int(covars[b].nunique()),
            "cramers_v": v,
            # A perfectly nested pair (every level of a maps to one level of b)
            # cannot be separated statistically at all.
            "fully_nested_a_in_b": bool(
                (pd.crosstab(covars[a], covars[b]) > 0).sum(axis=1).max() == 1),
            "fully_nested_b_in_a": bool(
                (pd.crosstab(covars[b], covars[a]) > 0).sum(axis=1).max() == 1),
        })

    pairs_df = pd.DataFrame(pairs).sort_values("cramers_v", ascending=False)
    pairs_df["verdict"] = np.where(
        pairs_df["fully_nested_a_in_b"] | pairs_df["fully_nested_b_in_a"],
        "fully nested - effects NOT separable",
        np.where(pairs_df["cramers_v"] >= STRONG_V, "strongly entangled",
                 np.where(pairs_df["cramers_v"] >= 0.25, "moderately entangled", "weak")))

    man.write(matrix.reset_index().rename(columns={"index": "covariate"}),
              "covariate_association_matrix",
              "Square bias-corrected Cramer's V matrix between covariates",
              subdir=SUBDIR)
    man.write(pairs_df, "covariate_association_pairs",
              "Covariate pairs ranked by Cramer's V, with nesting and a verdict",
              subdir=SUBDIR)
    warn = pairs_df[pairs_df["verdict"] != "weak"]
    man.write(warn, "confound_warnings",
              "Covariate pairs entangled enough to confound downstream comparisons",
              subdir=SUBDIR)
    if not warn.empty:
        top = warn.iloc[0]
        log(f"  strongest entanglement: {top['covariate_a']} x {top['covariate_b']} "
            f"(V={top['cramers_v']:.3f}, {top['verdict']})")

    # -- how much QC variation does each covariate explain? ------------------
    if not qc.empty:
        rows = []
        for role, s in covars.items():
            for metric in qc.columns:
                rows.append({"covariate": role, "qc_metric": metric,
                             "eta_squared": eta_squared(qc[metric], s)})
        eta = pd.DataFrame(rows).sort_values("eta_squared", ascending=False)
        man.write(eta, "qc_variance_explained",
                  "eta-squared: fraction of each QC metric's variance explained "
                  "by each covariate (high = technical structure)", subdir=SUBDIR)

    for a, b in KEY_PAIRS:
        if a in covars and b in covars:
            tab = pd.crosstab(covars[a], covars[b])
            tab.index.name = a
            man.write(tab.reset_index(), f"crosstab_{a}_x_{b}",
                      f"Cell counts {a} x {b} -- read zero cells as unsupported contrasts",
                      subdir=SUBDIR)
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    for key, chem, ns in cli.dataset_variants(args):
        run(key, args, chem, ns)


if __name__ == "__main__":
    main()
