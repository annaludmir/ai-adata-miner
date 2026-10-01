#!/usr/bin/env python3
"""00 - Ground-truth inventory of each h5ad.

Reads the real file and writes what is *actually* there: obs/var columns with
dtypes and cardinality, obsm/varm/layers/uns keys, and matrix shape.  Run this
first -- every later script trusts the file, not schemas/*.json, and this is
the step that proves the two agree.

Outputs
  obs_columns.csv      one row per .obs column: dtype, n_unique, missingness
  var_columns.csv      one row per .var column
  element_keys.csv     obsm / varm / layers / uns / obsp inventory
  dataset_summary.csv  shape, encoding, and schema-agreement check
  <dataset>_schema_observed.json  regenerated schema next to the stored one
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

import config
from lib import cli
from lib.io_utils import (Manifest, dump_json, list_h5ad_keys, log, read_obs,
                          read_var)

SCRIPT = "00_inventory_and_schema"


def describe_frame(df: pd.DataFrame, kind: str) -> pd.DataFrame:
    rows = []
    for col in df.columns:
        s = df[col]
        n_unique = int(s.nunique(dropna=True))
        # Preview values without dumping a 600k-category list into a CSV cell.
        try:
            sample = "|".join(map(str, pd.unique(s.dropna())[:8]))[:300]
        except Exception:
            sample = ""
        rows.append({
            "kind": kind,
            "column": col,
            "dtype": str(s.dtype),
            "n_unique": n_unique,
            "n_missing": int(s.isna().sum()),
            "frac_missing": float(s.isna().mean()),
            "is_categorical": isinstance(s.dtype, pd.CategoricalDtype),
            "is_numeric": bool(pd.api.types.is_numeric_dtype(s)),
            "constant": n_unique <= 1,
            "sample_values": sample,
        })
    return pd.DataFrame(rows)


def run(key: str, args, chem: str | None = None,
        ns: str | None = None) -> None:
    ns = ns or key
    cli.banner(SCRIPT, key, chem)
    path = cli.resolve_h5ad(key)
    man = Manifest(ns, SCRIPT)
    meta = config.dataset(key)

    log("reading .obs (no X touched)...")
    obs = read_obs(path)
    log(f"  obs: {obs.shape[0]:,} cells x {obs.shape[1]} columns")
    var = read_var(path)
    log(f"  var: {var.shape[0]:,} genes x {var.shape[1]} columns")
    keys = list_h5ad_keys(path)

    man.write(describe_frame(obs, "obs"), "obs_columns",
              "One row per .obs column: dtype, cardinality, missingness, sample values",
              subdir="00_inventory")
    man.write(describe_frame(var, "var"), "var_columns",
              "One row per .var column", subdir="00_inventory")

    elements = pd.DataFrame(
        [{"group": g, "key": k} for g, ks in keys.items()
         if g in ("obsm", "varm", "layers", "uns", "obsp") for k in ks])
    if elements.empty:
        elements = pd.DataFrame(columns=["group", "key"])
    man.write(elements, "element_keys",
              "Inventory of obsm/varm/layers/uns/obsp keys available for analysis",
              subdir="00_inventory")

    shape = keys.get("shape") or [obs.shape[0], var.shape[0]]
    summary = pd.DataFrame([{
        "dataset": key,
        "label": meta["label"],
        "h5ad": str(path),
        "n_obs": int(shape[0]),
        "n_vars": int(shape[1]),
        "X_encoding": keys["X_encoding"][0],
        "n_obs_columns": obs.shape[1],
        "n_var_columns": var.shape[1],
        "var_index_name": str(var.index.name),
        "gene_id_space": meta["gene_id_space"],
        "obsm_keys": "|".join(keys["obsm"]),
        "varm_keys": "|".join(keys["varm"]),
        "layers": "|".join(keys["layers"]),
        "matches_expected_n_obs": int(shape[0]) == meta["expected_n_obs"],
        "matches_expected_n_vars": int(shape[1]) == meta["expected_n_vars"],
    }])
    man.write(summary, "dataset_summary",
              "Shape, encoding and agreement with the stored schema",
              subdir="00_inventory")

    if not bool(summary.loc[0, "matches_expected_n_obs"]):
        log(f"  WARNING: n_obs {shape[0]:,} != schema's {meta['expected_n_obs']:,}")
    if not bool(summary.loc[0, "matches_expected_n_vars"]):
        log(f"  WARNING: n_vars {shape[1]:,} != schema's {meta['expected_n_vars']:,}")

    observed = {
        "dataset": key, "n_obs": int(shape[0]), "n_vars": int(shape[1]),
        "obs_index": {"name": str(obs.index.name),
                      "sample_values": list(map(str, obs.index[:5]))},
        "var_index": {"name": str(var.index.name),
                      "sample_values": list(map(str, var.index[:5]))},
        "obs": {c: {"dtype": str(obs[c].dtype), "n_unique": int(obs[c].nunique())}
                for c in obs.columns},
        "var": {c: {"dtype": str(var[c].dtype), "n_unique": int(var[c].nunique())}
                for c in var.columns},
        **{f"{g}_keys": keys[g] for g in ("obsm", "varm", "uns")},
        "layers": keys["layers"],
    }
    dump_json(observed, config.SCHEMA_DIR / f"{key}_schema_observed.json")
    log(f"  regenerated schema -> schemas/{key}_schema_observed.json")
    man.flush()


def main() -> None:
    args = cli.build_parser(__doc__).parse_args()
    if getattr(args, "chemistry", "all") != "all":
        log("  (this step describes the file as a whole; chemistry split not applied)")
    for key in cli.selected_datasets(args):
        run(key, args, None, key)


if __name__ == "__main__":
    main()
