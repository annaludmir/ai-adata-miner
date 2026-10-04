"""Thin, defensive h5ad access + CSV export bookkeeping.

Design notes
------------
* Nothing opens a full AnnData.  `read_obs` pulls only the /obs group, and
  `XReader` slices CSR row-blocks straight out of HDF5, so streaming the
  1.67M-cell matrix never materialises obsm (Factors alone would be ~GB).
* Column names are never hard-coded; `resolve_role` maps a semantic role to
  whichever candidate column the file actually has (see config.COLUMN_ROLES).
* Every CSV written goes through `Manifest`, which drops a machine-readable
  index at csv_exports/<dataset>/_manifest.csv.  Step 3 reads that manifest to
  discover what is available instead of guessing filenames.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import config

__all__ = [
    "log", "read_obs", "read_var", "read_elem_at", "list_h5ad_keys", "XReader",
    "resolve_role", "resolve_cluster_columns", "resolve_qc_frame",
    "coerce_numeric", "gene_frame", "harmonise_labels", "Manifest", "chemistry_mask", "load_obs",
    "add_derived_obs_columns",
]


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# h5ad reading
# ---------------------------------------------------------------------------
def _read_elem():
    """anndata moved read_elem twice; support every layout we might meet."""
    try:
        from anndata.io import read_elem            # anndata >= 0.11
    except ImportError:
        try:
            from anndata.experimental import read_elem   # 0.8 - 0.10
        except ImportError:
            from anndata._io.specs import read_elem      # fallback
    return read_elem


def _h5py():
    import h5py
    return h5py


def read_obs(path: Path | str) -> pd.DataFrame:
    """Read .obs only -- no X, no obsm. Seconds even on the 1.67M-cell file."""
    h5py, read_elem = _h5py(), _read_elem()
    with h5py.File(str(path), "r") as f:
        obs = read_elem(f["obs"])
    return pd.DataFrame(obs)


def read_var(path: Path | str) -> pd.DataFrame:
    h5py, read_elem = _h5py(), _read_elem()
    with h5py.File(str(path), "r") as f:
        var = read_elem(f["var"])
    return pd.DataFrame(var)


def read_elem_at(path: Path | str, key: str):
    """Read one arbitrary element, e.g. 'obsm/Factors' or 'varm/Loadings'."""
    h5py, read_elem = _h5py(), _read_elem()
    with h5py.File(str(path), "r") as f:
        if key not in f:
            return None
        return read_elem(f[key])


def list_h5ad_keys(path: Path | str) -> dict[str, list[str]]:
    """Inventory obsm/varm/layers/uns/obsp without loading anything."""
    h5py = _h5py()
    out: dict[str, list[str]] = {}
    with h5py.File(str(path), "r") as f:
        for grp in ("obsm", "varm", "layers", "uns", "obsp"):
            out[grp] = sorted(f[grp].keys()) if grp in f else []
        x = f.get("X")
        out["X_encoding"] = [x.attrs.get("encoding-type", "dense") if x is not None else "missing"]
        out["shape"] = list(_shape_of(f, "X")) if x is not None else []
    return out


def _shape_of(f, key: str) -> tuple[int, int]:
    node = f[key]
    if "shape" in node.attrs:
        return tuple(int(v) for v in node.attrs["shape"])
    return tuple(int(v) for v in node.shape)


class XReader:
    """Row-chunk reader for /X or /layers/<name>, sparse-CSR or dense.

    Used as a context manager:

        with XReader(path) as xr:
            for start, stop, chunk in xr.iter_chunks(50_000):
                ...
    """

    def __init__(self, path: Path | str, key: str = "X"):
        self.path = str(path)
        self.key = key
        self._f = None

    def __enter__(self) -> "XReader":
        self._f = _h5py().File(self.path, "r")
        if self.key not in self._f:
            self._f.close()
            raise KeyError(f"{self.key!r} not present in {self.path}")
        node = self._f[self.key]
        self.encoding = node.attrs.get("encoding-type", None)
        if self.encoding is None:
            self.encoding = "csr_matrix" if isinstance(node, _h5py().Group) else "array"
        if isinstance(self.encoding, bytes):
            self.encoding = self.encoding.decode()
        self.shape = _shape_of(self._f, self.key)
        if self.encoding == "csc_matrix":
            raise NotImplementedError(
                f"{self.key} is CSC; row-chunked streaming needs CSR. "
                "Re-save the file with a CSR X, or load it with anndata in backed mode.")
        return self

    def __exit__(self, *exc) -> None:
        if self._f is not None:
            self._f.close()
            self._f = None

    @property
    def n_obs(self) -> int:
        return int(self.shape[0])

    @property
    def n_vars(self) -> int:
        return int(self.shape[1])

    def read_rows(self, start: int, stop: int) -> sp.csr_matrix:
        node = self._f[self.key]
        stop = min(stop, self.shape[0])
        if self.encoding == "array":
            return sp.csr_matrix(np.asarray(node[start:stop, :], dtype=np.float64))
        indptr = node["indptr"][start:stop + 1]
        lo, hi = int(indptr[0]), int(indptr[-1])
        data = np.asarray(node["data"][lo:hi], dtype=np.float64)
        indices = np.asarray(node["indices"][lo:hi], dtype=np.int64)
        return sp.csr_matrix(
            (data, indices, indptr - lo),
            shape=(stop - start, self.shape[1]),
        )

    def iter_chunks(self, chunk_size: int = config.CHUNK_SIZE):
        for start in range(0, self.shape[0], chunk_size):
            stop = min(start + chunk_size, self.shape[0])
            yield start, stop, self.read_rows(start, stop)


# ---------------------------------------------------------------------------
# Column resolution
# ---------------------------------------------------------------------------
def resolve_role(obs: pd.DataFrame, role: str) -> str | None:
    """First obs column matching a semantic role, or None if the file lacks it."""
    for cand in config.COLUMN_ROLES.get(role, []):
        if cand in obs.columns:
            return cand
    return None


def resolve_cluster_columns(obs: pd.DataFrame) -> list[str]:
    """Every clustering column present -- alternative partitions are all profiled."""
    return [c for c in config.CLUSTER_ROLES if c in obs.columns]


def coerce_numeric(s: pd.Series) -> pd.Series:
    """Force a column to float.

    human_dev stores total_genes / total_UMIs as *categorical strings*; a naive
    .mean() on those silently fails or sorts lexically. Empty strings (seen in
    cortex's Cellconc / Targetnumcells) become NaN rather than 0.
    """
    if isinstance(s.dtype, pd.CategoricalDtype):
        s = s.astype("object")
    if s.dtype == object or pd.api.types.is_string_dtype(s):
        s = s.replace({"": None, "N/A": None, "NA": None, "nan": None, "None": None})
    if s.dtype == bool:
        return s.astype(np.float64)
    # Always land on float64: a mix of int64 / Int64 / float across datasets
    # otherwise bites downstream wherever NaN has to be representable.
    return pd.to_numeric(s, errors="coerce").astype(np.float64)


def resolve_qc_frame(obs: pd.DataFrame) -> pd.DataFrame:
    """Numeric QC metrics present in this file, renamed to canonical names."""
    cols: dict[str, pd.Series] = {}
    for canonical, candidates in config.QC_NUMERIC_ROLES.items():
        for cand in candidates:
            if cand in obs.columns:
                coerced = coerce_numeric(obs[cand])
                if coerced.notna().any():
                    cols[canonical] = coerced
                break
    return pd.DataFrame(cols, index=obs.index)


def harmonise_labels(s: pd.Series, role: str) -> pd.Series:
    """Apply cross-dataset label synonyms (e.g. 'PostM' -> 'Post-M')."""
    mapping = config.LABEL_SYNONYMS.get(role)
    if not mapping:
        return s.astype("object")
    return s.astype("object").map(lambda v: mapping.get(v, v))


def add_derived_obs_columns(obs: pd.DataFrame, dataset_key: str) -> pd.DataFrame:
    """Add analysis-ready derived columns without mutating the caller's frame.

    * age_pcw      -- Age parsed to float and rounded (human_dev stores it as a
                      categorical string like '8.100000381469727').
    * age_bin_week -- integer post-conception week, for coarse trajectories.
    * cyclephase_h -- cell-cycle phase harmonised across the two datasets.
    """
    out = obs.copy()
    age_col = resolve_role(obs, "age")
    if age_col is not None:
        age = coerce_numeric(obs[age_col])
        out["age_pcw"] = age.round(2)
        out["age_bin_week"] = np.floor(age).astype("Int64")
    phase_col = resolve_role(obs, "cyclephase")
    if phase_col is not None:
        out["cyclephase_h"] = harmonise_labels(obs[phase_col], "cyclephase")
    out["dataset"] = dataset_key
    return out


def normalise_donor(d) -> str:
    """cortex writes 'XHU:1966:307' where human_dev writes 'XHU:307'."""
    import re
    return re.sub(r"^(XHU|XDD):\d+:(\d+)$", r"\1:\2", str(d))


def load_exclusions(dataset_key: str) -> pd.DataFrame:
    """Exclusion rules for one dataset ('*' in the dataset column = every dataset)."""
    cols = ["dataset", "role", "value", "reason"]
    src = str(config.EXCLUSIONS_FILE or "")
    if src.lower() in ("", "none") or not Path(src).exists():
        if src.lower() not in ("", "none"):
            log(f"  exclusions file {src} not found -- no cells excluded")
        return pd.DataFrame(columns=cols)
    rules = pd.read_csv(src, dtype=str, comment="#").fillna("")
    missing = set(cols) - set(rules.columns)
    if missing:
        raise ValueError(f"{src} is missing columns {sorted(missing)}")
    return rules[rules["dataset"].isin([dataset_key, "*"])].reset_index(drop=True)


def exclusion_mask(obs: pd.DataFrame, dataset_key: str,
                   within: np.ndarray | None = None) -> tuple[np.ndarray, pd.DataFrame]:
    """(keep, report): which cells survive the exclusion rules, and what each removed.

    A rule's role is a semantic role from config.COLUMN_ROLES (age, donor,
    sample, region, ...) or a raw obs column name. Ages match numerically
    (to 0.01 pcw, so '5' matches the categorical string '5.0'); donors match
    after normalising cortex-style IDs; everything else matches as a string.
    `within` (e.g. the chemistry mask) restricts the counts in the report and
    the log to the cells actually being analysed.
    """
    rules = load_exclusions(dataset_key)
    keep = np.ones(len(obs), dtype=bool)
    scope = np.ones(len(obs), dtype=bool) if within is None else np.asarray(within, bool)
    report = []
    for r in rules.itertuples():
        col = resolve_role(obs, r.role) or (r.role if r.role in obs.columns else None)
        if col is None:
            log(f"  WARNING: exclusion rule {r.role}={r.value} -- no such column in "
                f"{dataset_key}; rule ignored")
            hit = np.zeros(len(obs), dtype=bool)
        elif r.role == "age":
            age = coerce_numeric(obs[col]).to_numpy()
            hit = np.abs(age - float(r.value)) < 0.01
        elif r.role == "donor":
            hit = (obs[col].astype(str).map(normalise_donor)
                   == normalise_donor(r.value)).to_numpy()
        else:
            hit = (obs[col].astype(str) == str(r.value)).to_numpy()
        report.append({"dataset": dataset_key, "role": r.role, "column": col or "",
                       "value": r.value, "reason": r.reason,
                       "n_cells_matched": int((hit & scope).sum()),
                       "n_cells_removed": int((hit & keep & scope).sum())})
        keep &= ~hit
    if report and (~keep & scope).any():
        log(f"  exclusions ({config.EXCLUSIONS_FILE}): removed {int((~keep & scope).sum()):,} cells -- "
            + "; ".join(f"{x['role']}={x['value']}: {x['n_cells_removed']:,}" for x in report))
    return keep, pd.DataFrame(report, columns=["dataset", "role", "column", "value", "reason",
                                               "n_cells_matched", "n_cells_removed"])


def load_obs(path, dataset_key: str, chemistry: str | None = None,
             limit_cells: int | None = None):
    """Read .obs, apply the chemistry filter and exclusions, add derived columns.

    Returns (obs, keep) where `keep` is a boolean array over the rows that were
    read, so callers can subset a positionally-aligned obsm matrix the same way.
    Returns (None, None) when a stratified run was requested but the file has no
    chemistry column -- the caller must skip rather than analyse every cell and
    mislabel it as one chemistry.
    """
    obs = read_obs(path)
    if limit_cells:
        obs = obs.iloc[:limit_cells]
    mask = chemistry_mask(obs, chemistry)
    if mask is None:
        return None, None
    excl_keep, report = exclusion_mask(obs, dataset_key, within=mask.to_numpy())
    keep = mask.to_numpy() & excl_keep
    out = add_derived_obs_columns(obs.loc[keep], dataset_key)
    out.attrs["exclusions"] = report     # what was removed, for the record
    return out, keep


def load_group_matrix(path, min_cells: int | None = None) -> pd.DataFrame:
    """Read a 09_pseudobulk gene x group table, keeping groups with enough cells.

    Uses the sibling `<grouping>__group_summary.csv`. Groups below `min_cells`
    (default config.MIN_CELLS_PER_GROUP) are dropped -- this also removes the
    all-zero columns that exports from before the empty-group fix contain.
    """
    path = Path(path)
    df = pd.read_csv(path)
    df = df.set_index(df.columns[0])
    df.index.name = "gene"
    df = df.apply(pd.to_numeric, errors="coerce")
    min_cells = config.MIN_CELLS_PER_GROUP if min_cells is None else min_cells
    grouping = path.name.split("__")[0]
    spath = path.with_name(f"{grouping}__group_summary.csv")
    if spath.exists():
        gs = pd.read_csv(spath)
        good = set(gs.loc[gs["n_cells"] >= min_cells, "group"].astype(str))
        keep = [c for c in df.columns if str(c) in good]
    else:
        keep = [c for c in df.columns if df[c].fillna(0).abs().sum() > 0]
    dropped = df.shape[1] - len(keep)
    if dropped:
        log(f"  {grouping}: dropped {dropped}/{df.shape[1]} groups with "
            f"< {min_cells} cells")
    return df[keep]


def chemistry_mask(obs: pd.DataFrame, chemistry: str | None) -> pd.Series | None:
    """Boolean mask selecting one chemistry's cells.

    Returns an all-True mask when `chemistry` is None (pooled), and None when a
    stratified run was asked for but the file has no chemistry column -- the
    caller must then skip rather than silently analyse everything and label it
    'v2'.
    """
    if chemistry is None:
        return pd.Series(True, index=obs.index)
    col = resolve_role(obs, "chemistry")
    if col is None:
        return None
    return obs[col].astype(str).str.strip() == chemistry


def gene_frame(var: pd.DataFrame, dataset_key: str) -> pd.DataFrame:
    """Unified gene table: symbol, accession, and a version-stripped accession.

    The two files live in different id spaces -- cortex is indexed by symbol
    with unversioned accessions, human_dev by *versioned* accession.  Stripping
    the version gives the one key that joins them.
    """
    g = pd.DataFrame(index=var.index.copy())
    g.index.name = "var_index"

    symbol_col = next((c for c in config.VAR_SYMBOL_CANDIDATES if c in var.columns), None)
    acc_col = next((c for c in config.VAR_ACCESSION_CANDIDATES if c in var.columns), None)

    # Guard against human_dev's mislabelled `gene_symbol` column, which holds
    # accessions. If the chosen symbol column looks like ENSG ids, reject it.
    def _looks_like_accession(col: str) -> bool:
        return var[col].astype(str).head(50).str.startswith("ENSG").mean() > 0.5

    if symbol_col is not None and _looks_like_accession(symbol_col):
        symbol_col = next(
            (c for c in config.VAR_SYMBOL_CANDIDATES
             if c in var.columns and c != symbol_col and not _looks_like_accession(c)),
            None)

    idx_str = var.index.astype(str)
    g["symbol"] = var[symbol_col].astype(str) if symbol_col else idx_str
    g["accession"] = var[acc_col].astype(str) if acc_col else (
        idx_str if idx_str.str.startswith("ENSG").mean() > 0.5 else pd.NA)
    g["accession_base"] = (
        g["accession"].astype(str).str.split(".").str[0].replace("nan", pd.NA)
        if g["accession"].notna().any() else pd.NA)
    for extra in ("Chromosome", "Start", "End", "GeneTotalUMIs"):
        if extra in var.columns:
            g[extra] = var[extra].values
    g["dataset"] = dataset_key
    return g


# ---------------------------------------------------------------------------
# CSV export bookkeeping
# ---------------------------------------------------------------------------
class Manifest:
    """Records every CSV written so step 3 can discover outputs programmatically."""

    def __init__(self, dataset_key: str, script: str):
        self.dataset_key = dataset_key
        self.script = script
        self.rows: list[dict] = []
        self.root = config.out_dir(dataset_key)

    def write(self, df: pd.DataFrame, name: str, description: str,
              subdir: str | None = None, index: bool = False,
              float_format: str | None = "%.6g") -> Path:
        target_dir = config.out_dir(self.dataset_key, subdir)
        path = target_dir / (name if name.endswith(".csv") else f"{name}.csv")
        df.to_csv(path, index=index, float_format=float_format)
        self.rows.append({
            "csv": str(path.relative_to(self.root.parent)),
            "dataset": self.dataset_key,
            "script": self.script,
            "n_rows": int(df.shape[0]),
            "n_cols": int(df.shape[1]),
            "columns": "|".join(map(str, df.columns[:40])),
            "index_name": str(df.index.name) if index else "",
            "description": description,
            "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        })
        log(f"  wrote {path.name}  ({df.shape[0]} x {df.shape[1]})  {description}")
        return path

    def flush(self) -> None:
        """Merge this run's entries into csv_exports/<dataset>/_manifest.csv."""
        if not self.rows:
            log("  no CSVs written")
            return
        path = self.root / "_manifest.csv"
        new = pd.DataFrame(self.rows)
        if path.exists():
            old = pd.read_csv(path)
            old = old[~old["csv"].isin(new["csv"])]
            new = pd.concat([old, new], ignore_index=True)
        new.sort_values(["script", "csv"]).to_csv(path, index=False)
        log(f"  manifest updated: {path} ({len(new)} entries)")


def require_file(path: Path) -> Path:
    if not Path(path).exists():
        raise SystemExit(
            f"Input file not found: {path}\n"
            f"Set AI_ADATA_DATA_ROOT if the data lives elsewhere.")
    return Path(path)


def dump_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))
