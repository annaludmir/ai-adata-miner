"""Streaming group aggregation over an expression matrix.

The whole-brain file is 1.67M cells x 59k genes, so nothing here ever
materialises the matrix.  Cells are grouped by an integer code vector and
accumulated chunk by chunk via a sparse one-hot product:

    group_sums += onehot(codes).T @ X_chunk        # (n_groups x n_genes)

One pass yields three complementary views of every group:
  * summed raw counts  -> true pseudobulk, the correct input for edgeR/DESeq2
  * mean log1p(CP10K)  -> comparable across groups of different depth
  * detection fraction -> how many cells express the gene at all

Depends only on numpy/scipy/pandas so it can be tested without an h5ad.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

__all__ = ["group_codes", "onehot_csr", "GroupAggregator", "contingency",
           "combine_keys"]

MISSING_CODE = -1


def combine_keys(df: pd.DataFrame, columns: list[str], sep: str = " | ") -> pd.Series:
    """Join several obs columns into one composite group label."""
    parts = [df[c].astype("string").fillna("NA") for c in columns]
    out = parts[0]
    for p in parts[1:]:
        out = out + sep + p
    return out.astype("string")


def group_codes(values: pd.Series) -> tuple[np.ndarray, list[str]]:
    """Map a label series to contiguous integer codes.

    Missing labels get MISSING_CODE (-1) and are excluded from every
    accumulator rather than silently folded into a bogus group.
    """
    cat = pd.Categorical(values.astype("object"))
    codes = np.asarray(cat.codes, dtype=np.int64)
    return codes, [str(c) for c in cat.categories]


def drop_empty_levels(codes: np.ndarray, levels: list[str]
                      ) -> tuple[np.ndarray, list[str]]:
    """Re-code so only levels with at least one kept cell remain.

    Levels come from the whole file, but a chemistry stratum masks the other
    chemistry's cells to -1; without this, every age, donor or sample that
    exists only in the other chemistry is exported as an all-zero column.
    """
    codes = np.asarray(codes, dtype=np.int64)
    present = np.unique(codes[codes >= 0])
    if len(present) == len(levels):
        return codes, levels
    remap = np.full(len(levels), -1, dtype=np.int64)
    remap[present] = np.arange(len(present))
    new = np.where(codes >= 0, remap[np.clip(codes, 0, None)], -1)
    return new, [levels[i] for i in present]


def onehot_csr(codes: np.ndarray, n_groups: int) -> sp.csr_matrix:
    """(n_cells x n_groups) indicator matrix; rows with code -1 are all-zero."""
    codes = np.asarray(codes, dtype=np.int64)
    keep = codes >= MISSING_CODE + 1
    rows = np.nonzero(keep)[0]
    cols = codes[keep]
    data = np.ones(rows.size, dtype=np.float64)
    return sp.csr_matrix((data, (rows, cols)), shape=(codes.size, n_groups))


def _row_normalise(X: sp.csr_matrix, totals: np.ndarray, target_sum: float):
    """Scale each row to `target_sum` counts, then log1p. Zero-count rows stay zero."""
    scale = np.zeros_like(totals, dtype=np.float64)
    nz = totals > 0
    scale[nz] = target_sum / totals[nz]
    out = sp.diags(scale) @ X
    out = out.tocsr()
    out.data = np.log1p(out.data)
    return out


class GroupAggregator:
    """Accumulate per-group expression statistics across chunks of cells."""

    def __init__(self, n_groups: int, n_genes: int, target_sum: float = 1e4):
        self.n_groups = int(n_groups)
        self.n_genes = int(n_genes)
        self.target_sum = float(target_sum)
        self.sum_counts = np.zeros((self.n_groups, self.n_genes), dtype=np.float64)
        self.sum_lognorm = np.zeros((self.n_groups, self.n_genes), dtype=np.float64)
        self.n_detected = np.zeros((self.n_groups, self.n_genes), dtype=np.float64)
        self.n_cells = np.zeros(self.n_groups, dtype=np.float64)
        self.sum_umis = np.zeros(self.n_groups, dtype=np.float64)
        self.n_cells_skipped = 0

    def update(self, X_chunk, codes_chunk: np.ndarray,
               totals_chunk: np.ndarray | None = None,
               gene_mask: np.ndarray | None = None) -> None:
        """Add one chunk of cells.

        `totals_chunk` is the per-cell UMI total over *all* genes and must be
        supplied when `gene_mask` subsets the matrix, otherwise CP10K would be
        normalised against a truncated library size.
        """
        X = sp.csr_matrix(X_chunk) if not sp.isspmatrix_csr(X_chunk) else X_chunk
        codes = np.asarray(codes_chunk, dtype=np.int64)
        if codes.size != X.shape[0]:
            raise ValueError(f"codes ({codes.size}) != rows in chunk ({X.shape[0]})")

        if totals_chunk is None:
            totals = np.asarray(X.sum(axis=1)).ravel()
        else:
            totals = np.asarray(totals_chunk, dtype=np.float64).ravel()

        if gene_mask is not None:
            X = X[:, gene_mask]
        if X.shape[1] != self.n_genes:
            raise ValueError(f"chunk has {X.shape[1]} genes, expected {self.n_genes}")

        self.n_cells_skipped += int((codes == MISSING_CODE).sum())
        G = onehot_csr(codes, self.n_groups)      # (cells x groups)
        Gt = G.T.tocsr()

        self.sum_counts += np.asarray((Gt @ X).todense())
        self.sum_lognorm += np.asarray((Gt @ _row_normalise(X, totals, self.target_sum)).todense())

        binarised = X.copy()
        binarised.data = np.ones_like(binarised.data, dtype=np.float64)
        self.n_detected += np.asarray((Gt @ binarised).todense())

        self.n_cells += np.asarray(G.sum(axis=0)).ravel()
        self.sum_umis += np.asarray(Gt @ totals).ravel()

    # -- outputs -----------------------------------------------------------
    def _frame(self, matrix: np.ndarray, groups: list[str], genes: list[str]) -> pd.DataFrame:
        """Genes as rows, groups as columns: the orientation CSV readers want."""
        return pd.DataFrame(matrix.T, index=pd.Index(genes, name="gene"), columns=groups)

    def pseudobulk_counts(self, groups, genes) -> pd.DataFrame:
        return self._frame(self.sum_counts, groups, genes)

    def mean_lognorm(self, groups, genes) -> pd.DataFrame:
        denom = np.where(self.n_cells[:, None] > 0, self.n_cells[:, None], np.nan)
        return self._frame(self.sum_lognorm / denom, groups, genes)

    def detection_fraction(self, groups, genes) -> pd.DataFrame:
        denom = np.where(self.n_cells[:, None] > 0, self.n_cells[:, None], np.nan)
        return self._frame(self.n_detected / denom, groups, genes)

    def cpm(self, groups, genes) -> pd.DataFrame:
        """Counts per million of the pseudobulk library -- depth-independent."""
        lib = self.sum_counts.sum(axis=1, keepdims=True)
        lib = np.where(lib > 0, lib, np.nan)
        return self._frame(self.sum_counts / lib * 1e6, groups, genes)

    def group_summary(self, groups: list[str]) -> pd.DataFrame:
        with np.errstate(divide="ignore", invalid="ignore"):
            mean_umis = np.where(self.n_cells > 0, self.sum_umis / self.n_cells, np.nan)
            genes_det = np.where(self.n_cells > 0,
                                 (self.n_detected > 0).sum(axis=1), np.nan)
        return pd.DataFrame({
            "group": groups,
            "n_cells": self.n_cells.astype(np.int64),
            "total_umis": self.sum_umis,
            "mean_umis_per_cell": mean_umis,
            "n_genes_detected_in_group": genes_det,
        })


def contingency(df: pd.DataFrame, row_col: str, col_col: str) -> pd.DataFrame:
    """Cell counts cross-tabulated, with NaN rows dropped."""
    sub = df[[row_col, col_col]].dropna()
    return pd.crosstab(sub[row_col].astype(str), sub[col_col].astype(str))
