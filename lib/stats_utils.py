"""Pure numeric helpers: diversity, association, specificity.

Deliberately free of h5py/anndata/scanpy so the statistics can be unit-tested
anywhere.  Everything here takes and returns numpy/pandas objects.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "shannon_entropy", "normalized_entropy", "simpson_index", "gini",
    "cramers_v", "chi2_standardised_residuals", "log2_observed_expected",
    "tau_specificity", "entropy_specificity", "corr_matrix",
]


# ---------------------------------------------------------------------------
# Diversity of a composition vector
# ---------------------------------------------------------------------------
def _as_proportions(counts) -> np.ndarray:
    c = np.asarray(counts, dtype=float).ravel()
    c = c[np.isfinite(c)]
    c = c[c > 0]
    total = c.sum()
    return c / total if total > 0 else np.array([])


def shannon_entropy(counts, base: float = np.e) -> float:
    """Shannon entropy of a count vector. 0 for a pure / empty composition."""
    p = _as_proportions(counts)
    if p.size == 0:
        return float("nan")
    return float(-(p * np.log(p)).sum() / np.log(base))


def normalized_entropy(counts) -> float:
    """Shannon entropy divided by log(number of observed categories).

    Lands in [0, 1]: 0 = one category dominates, 1 = perfectly even.  Using the
    *observed* category count makes this comparable across groups that have
    different numbers of available categories.
    """
    p = _as_proportions(counts)
    if p.size <= 1:
        return 0.0 if p.size == 1 else float("nan")
    return float(-(p * np.log(p)).sum() / np.log(p.size))


def simpson_index(counts) -> float:
    """Probability two cells drawn at random come from different categories."""
    p = _as_proportions(counts)
    if p.size == 0:
        return float("nan")
    return float(1.0 - (p ** 2).sum())


def gini(values) -> float:
    """Gini coefficient of a non-negative vector (0 = even, ->1 = concentrated)."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0 or np.all(v <= 0):
        return float("nan")
    v = np.sort(np.clip(v, 0, None))
    n = v.size
    cum = np.cumsum(v)
    return float((n + 1 - 2 * (cum / cum[-1]).sum()) / n)


# ---------------------------------------------------------------------------
# Association between two categorical variables
# ---------------------------------------------------------------------------
def cramers_v(table: pd.DataFrame | np.ndarray, bias_correction: bool = True) -> float:
    """Cramer's V for a contingency table.

    With `bias_correction` this is the Bergsma-Wicher corrected V, which stops
    small or sparse tables from reporting spurious near-1 association.  Returns
    NaN for degenerate tables (a single row or column, or no observations).
    """
    obs = np.asarray(table, dtype=float)
    n = obs.sum()
    if n <= 0 or obs.shape[0] < 2 or obs.shape[1] < 2:
        return float("nan")

    row = obs.sum(axis=1, keepdims=True)
    col = obs.sum(axis=0, keepdims=True)
    expected = row @ col / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.where(expected > 0, (obs - expected) ** 2 / expected, 0.0).sum()

    phi2 = chi2 / n
    r, k = obs.shape
    if not bias_correction:
        denom = min(r - 1, k - 1)
        return float(np.sqrt(phi2 / denom)) if denom > 0 else float("nan")

    phi2_corr = max(0.0, phi2 - (k - 1) * (r - 1) / max(n - 1, 1))
    r_corr = r - (r - 1) ** 2 / max(n - 1, 1)
    k_corr = k - (k - 1) ** 2 / max(n - 1, 1)
    denom = min(r_corr - 1, k_corr - 1)
    if denom <= 0:
        return float("nan")
    return float(np.sqrt(phi2_corr / denom))


def chi2_standardised_residuals(table: pd.DataFrame) -> pd.DataFrame:
    """Standardised Pearson residuals: which cells of the table are enriched.

    Values are approximately z-scores under independence, so |r| > 2 flags a
    group/category combination that departs from what marginals predict.
    """
    obs = table.to_numpy(dtype=float)
    n = obs.sum()
    if n <= 0:
        return pd.DataFrame(np.nan, index=table.index, columns=table.columns)
    row = obs.sum(axis=1, keepdims=True)
    col = obs.sum(axis=0, keepdims=True)
    expected = row @ col / n
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = np.sqrt(expected * (1 - row / n) * (1 - col / n))
        res = np.where(denom > 0, (obs - expected) / denom, np.nan)
    return pd.DataFrame(res, index=table.index, columns=table.columns)


def log2_observed_expected(table: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """log2(observed / expected) per cell of a contingency table.

    More interpretable than residuals for composition work: +1 means a cell
    type is twice as abundant in that group as the marginals would predict.
    """
    obs = table.to_numpy(dtype=float)
    n = obs.sum()
    if n <= 0:
        return pd.DataFrame(np.nan, index=table.index, columns=table.columns)
    expected = obs.sum(axis=1, keepdims=True) @ obs.sum(axis=0, keepdims=True) / n
    ratio = (obs + pseudocount) / (expected + pseudocount)
    return pd.DataFrame(np.log2(ratio), index=table.index, columns=table.columns)


# ---------------------------------------------------------------------------
# How specific is a gene / factor to one group?
# ---------------------------------------------------------------------------
def tau_specificity(matrix: np.ndarray) -> np.ndarray:
    """Yanai's tau per row of a (features x groups) non-negative matrix.

    0 = uniformly expressed across groups, 1 = restricted to a single group.
    Rows that are entirely zero return NaN rather than a misleading 0 or 1.
    """
    x = np.asarray(matrix, dtype=float)
    if x.ndim == 1:
        x = x[None, :]
    x = np.clip(np.nan_to_num(x, nan=0.0), 0, None)
    n_groups = x.shape[1]
    if n_groups < 2:
        return np.full(x.shape[0], np.nan)
    peak = x.max(axis=1)
    out = np.full(x.shape[0], np.nan)
    ok = peak > 0
    if ok.any():
        scaled = x[ok] / peak[ok][:, None]
        out[ok] = (1.0 - scaled).sum(axis=1) / (n_groups - 1)
    return out


def entropy_specificity(matrix: np.ndarray) -> np.ndarray:
    """1 - normalised entropy of each row treated as a distribution over groups."""
    x = np.clip(np.nan_to_num(np.asarray(matrix, dtype=float), nan=0.0), 0, None)
    if x.ndim == 1:
        x = x[None, :]
    n_groups = x.shape[1]
    if n_groups < 2:
        return np.full(x.shape[0], np.nan)
    totals = x.sum(axis=1, keepdims=True)
    out = np.full(x.shape[0], np.nan)
    ok = totals.ravel() > 0
    if ok.any():
        p = x[ok] / totals[ok]
        with np.errstate(divide="ignore", invalid="ignore"):
            ent = -(np.where(p > 0, p * np.log(p), 0.0)).sum(axis=1)
        out[ok] = 1.0 - ent / np.log(n_groups)
    return out


def corr_matrix(df: pd.DataFrame, method: str = "pearson",
                min_variance: float = 1e-12) -> pd.DataFrame:
    """Correlation matrix with zero-variance columns dropped up front.

    pandas would emit a RuntimeWarning and NaN column for constant inputs;
    dropping them keeps the exported CSV honest about what was comparable.
    """
    numeric = df.select_dtypes(include=[np.number])
    keep = [c for c in numeric.columns if numeric[c].var(skipna=True) > min_variance]
    if len(keep) < 2:
        return pd.DataFrame(dtype=float)
    return numeric[keep].corr(method=method)


# ---------------------------------------------------------------------------
# Comparing two partitions of the same cells
# ---------------------------------------------------------------------------
def _comb2(x: np.ndarray) -> np.ndarray:
    return x * (x - 1.0) / 2.0


def adjusted_rand_index(table: pd.DataFrame | np.ndarray) -> float:
    """ARI from a contingency table of two labellings of the same cells.

    0 = agreement no better than chance, 1 = identical partitions.  Lets us ask
    whether `louvain`, `leiden_scVI` and `Clusters` are really telling different
    stories or just relabelling the same structure.
    """
    n_ij = np.asarray(table, dtype=float)
    n = n_ij.sum()
    if n < 2:
        return float("nan")
    a, b = n_ij.sum(axis=1), n_ij.sum(axis=0)
    index = _comb2(n_ij).sum()
    exp = _comb2(a).sum() * _comb2(b).sum() / _comb2(np.array([n]))[0]
    mx = 0.5 * (_comb2(a).sum() + _comb2(b).sum())
    return float((index - exp) / (mx - exp)) if mx != exp else float("nan")


def normalized_mutual_info(table: pd.DataFrame | np.ndarray) -> float:
    """NMI with arithmetic-mean normalisation, from a contingency table."""
    n_ij = np.asarray(table, dtype=float)
    n = n_ij.sum()
    if n <= 0:
        return float("nan")
    p_ij = n_ij / n
    p_i = p_ij.sum(axis=1, keepdims=True)
    p_j = p_ij.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        mi = np.where(p_ij > 0, p_ij * np.log(p_ij / (p_i @ p_j)), 0.0).sum()
        h_u = -np.where(p_i > 0, p_i * np.log(p_i), 0.0).sum()
        h_v = -np.where(p_j > 0, p_j * np.log(p_j), 0.0).sum()
    denom = (h_u + h_v) / 2.0
    if denom <= 0:
        return 1.0 if mi == 0 else float("nan")
    return float(mi / denom)


def eta_squared(values: pd.Series, groups: pd.Series) -> float:
    """Fraction of a numeric column's variance explained by a categorical one.

    The continuous counterpart to Cramer's V: answers "how much of the
    mitochondrial fraction is explained by which donor the cell came from?"
    """
    df = pd.DataFrame({"v": pd.to_numeric(values, errors="coerce"),
                       "g": groups.astype("object")}).dropna()
    if df.empty or df["g"].nunique() < 2:
        return float("nan")
    grand = df["v"].mean()
    ss_total = ((df["v"] - grand) ** 2).sum()
    if ss_total <= 0:
        return float("nan")
    means = df.groupby("g", observed=True)["v"].agg(["mean", "count"])
    ss_between = (means["count"] * (means["mean"] - grand) ** 2).sum()
    return float(ss_between / ss_total)


__all__ += ["adjusted_rand_index", "normalized_mutual_info", "eta_squared"]
