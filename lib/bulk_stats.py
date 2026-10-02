"""Bulk-style expression statistics applied to the pseudobulk matrices.

Implements the classical RNA-seq workflow pieces -- TMM normalisation, PCA,
hierarchical clustering, K-means expression patterns and GSEA -- directly on
numpy/scipy so they can be unit-tested and so the pipeline keeps no heavyweight
dependency. Each function documents the assumption it makes, because several of
them are inherited from bulk RNA-seq and do not transfer unchanged to 10x data
(see `gene_length_note`).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "tmm_normalization_factors", "ma_statistics", "pca",
    "hierarchical_linkage", "cut_tree_clusters", "kmeans_clusters",
    "cluster_quality", "zscore_rows", "gsea_enrichment_score",
    "gsea_test", "benjamini_hochberg", "gene_length_note",
]

gene_length_note = """\
RPKM/FPKM is NOT applied anywhere in this pipeline, deliberately.

Both datasets are 10x UMI data. RPKM exists to correct the length bias of
full-length bulk protocols, where a longer transcript yields proportionally
more fragments. 10x counts one UMI per captured molecule from the 3' end, so
that bias is absent and dividing by gene length would introduce an artefact
rather than remove one. Within-sample comparisons across genes use detection
fraction and CPM; between-sample comparisons use TMM or CP10K."""


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------
def tmm_normalization_factors(counts: pd.DataFrame, ref_column: str | None = None,
                              logratio_trim: float = 0.3, sum_trim: float = 0.05,
                              min_counts: float = 0.0) -> pd.DataFrame:
    """Trimmed Mean of M-values scaling factors (Robinson & Oshlack 2010).

    CPM assumes every library samples the same underlying population; when a
    handful of genes dominate one group, that assumption fails and *every* other
    gene appears depressed. TMM estimates the composition bias from the trimmed
    mean of log-ratios against a reference library and corrects for it.

    `counts` is genes x samples of raw counts. Returns one row per sample with
    the library size, the raw factor, and the factor normalised to geometric
    mean 1 (which is how edgeR applies it).
    """
    X = counts.to_numpy(dtype=float)
    X = np.where(np.isfinite(X), X, 0.0)
    lib = X.sum(axis=0)
    samples = list(counts.columns)
    if X.shape[1] < 2 or not np.any(lib > 0):
        return pd.DataFrame({"sample": samples, "library_size": lib,
                             "tmm_factor_raw": np.nan, "tmm_factor": np.nan,
                             "effective_library_size": np.nan,
                             "n_genes_used": 0, "is_reference": False})

    # Reference = library whose upper-quartile CPM is closest to the mean of all.
    with np.errstate(invalid="ignore", divide="ignore"):
        cpm = np.where(lib > 0, X / lib, np.nan)
    uq = np.nanquantile(np.where(cpm > 0, cpm, np.nan), 0.75, axis=0)
    ref_idx = (samples.index(ref_column) if ref_column in samples
               else int(np.nanargmin(np.abs(uq - np.nanmean(uq)))))

    yr, nr = X[:, ref_idx], lib[ref_idx]
    factors, used = np.full(len(samples), np.nan), np.zeros(len(samples), dtype=int)
    for k in range(len(samples)):
        yk, nk = X[:, k], lib[k]
        if nk <= 0 or nr <= 0:
            continue
        keep = (yk > min_counts) & (yr > min_counts)
        if keep.sum() < 10:
            factors[k] = 1.0
            continue
        a_k, a_r = yk[keep] / nk, yr[keep] / nr
        M = np.log2(a_k / a_r)
        A = 0.5 * np.log2(a_k * a_r)
        # Inverse-variance weights from the delta method on the binomial.
        w = 1.0 / ((nk - yk[keep]) / (nk * yk[keep]) + (nr - yr[keep]) / (nr * yr[keep]))
        ok = np.isfinite(M) & np.isfinite(A) & np.isfinite(w)
        M, A, w = M[ok], A[ok], w[ok]
        if M.size < 10:
            factors[k] = 1.0
            continue
        # Double trim: on the log-ratio and on the average abundance.
        mlo, mhi = np.quantile(M, [logratio_trim, 1 - logratio_trim])
        alo, ahi = np.quantile(A, [sum_trim, 1 - sum_trim])
        sel = (M >= mlo) & (M <= mhi) & (A >= alo) & (A <= ahi)
        if sel.sum() == 0 or w[sel].sum() <= 0:
            factors[k] = 1.0
            continue
        factors[k] = 2.0 ** (np.sum(w[sel] * M[sel]) / np.sum(w[sel]))
        used[k] = int(sel.sum())

    finite = np.isfinite(factors) & (factors > 0)
    norm = factors.copy()
    if finite.any():
        norm[finite] = factors[finite] / np.exp(np.mean(np.log(factors[finite])))
    return pd.DataFrame({
        "sample": samples,
        "library_size": lib,
        "tmm_factor_raw": factors,
        "tmm_factor": norm,
        "effective_library_size": lib * norm,
        "n_genes_used": used,
        "is_reference": [i == ref_idx for i in range(len(samples))],
    })


def ma_statistics(counts: pd.DataFrame, ref_column: str | None = None) -> pd.DataFrame:
    """Per-sample MA-plot summary against a reference library.

    M is the log2 ratio and A the mean abundance. Under the normalisation
    assumption (no global shift) median M should sit at 0; a median far from
    zero is the bias the PDF's MA plots are drawn to expose, and the signal that
    CPM alone is not enough.
    """
    X = counts.to_numpy(dtype=float)
    lib = X.sum(axis=0)
    samples = list(counts.columns)
    if X.shape[1] < 2:
        return pd.DataFrame(columns=["sample", "median_M"])
    with np.errstate(invalid="ignore", divide="ignore"):
        uq = np.nanquantile(np.where(X > 0, X / np.where(lib > 0, lib, np.nan), np.nan),
                            0.75, axis=0)
    ref_idx = (samples.index(ref_column) if ref_column in samples
               else int(np.nanargmin(np.abs(uq - np.nanmean(uq)))))
    yr, nr = X[:, ref_idx], lib[ref_idx]
    rows = []
    for k, name in enumerate(samples):
        yk, nk = X[:, k], lib[k]
        keep = (yk > 0) & (yr > 0) & (nk > 0) & (nr > 0)
        if keep.sum() < 10:
            rows.append({"sample": name, "n_genes_compared": int(keep.sum())})
            continue
        M = np.log2((yk[keep] / nk) / (yr[keep] / nr))
        A = 0.5 * np.log2((yk[keep] / nk) * (yr[keep] / nr))
        rows.append({
            "sample": name, "reference_sample": samples[ref_idx],
            "n_genes_compared": int(keep.sum()),
            "median_M": float(np.median(M)), "mean_M": float(np.mean(M)),
            "iqr_M": float(np.subtract(*np.quantile(M, [0.75, 0.25]))),
            "frac_abs_M_gt1": float(np.mean(np.abs(M) > 1)),
            "median_A": float(np.median(A)),
            "bias_flag": bool(abs(float(np.median(M))) > 0.2),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Global structure
# ---------------------------------------------------------------------------
def pca(matrix: pd.DataFrame, n_components: int = 10, scale: bool = False):
    """PCA of samples (columns) over features (rows), via SVD.

    Returns (scores, loadings, variance table). Scores are samples x PCs, so
    this is the sample-relationship plot of the PDF's step 2; loadings say which
    genes drive each component.
    """
    X = matrix.to_numpy(dtype=float).T            # samples x genes
    X = np.nan_to_num(X, nan=0.0)
    keep = X.std(axis=0) > 0
    X, genes = X[:, keep], matrix.index.to_numpy()[keep]
    if X.shape[0] < 2 or X.shape[1] < 2:
        empty = pd.DataFrame()
        return empty, empty, empty
    X = X - X.mean(axis=0, keepdims=True)
    if scale:
        X = X / X.std(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    k = int(min(n_components, S.size))
    var = S ** 2 / max(X.shape[0] - 1, 1)
    ratio = var / var.sum() if var.sum() > 0 else np.full_like(var, np.nan)
    names = [f"PC{i+1}" for i in range(k)]
    scores = pd.DataFrame((U[:, :k] * S[:k]), index=matrix.columns, columns=names)
    scores.index.name = "sample"
    loadings = pd.DataFrame(Vt[:k].T, index=pd.Index(genes, name="gene"), columns=names)
    variance = pd.DataFrame({
        "component": names,
        "eigenvalue": var[:k],
        "proportion_of_variance": ratio[:k],
        "cumulative_proportion": np.cumsum(ratio)[:k],
    })
    return scores, loadings, variance


def hierarchical_linkage(matrix: pd.DataFrame, method: str = "average",
                         metric: str = "correlation"):
    """Hierarchical clustering of the columns. Returns (linkage, labels).

    Default average linkage on 1-Pearson distance: the PDF's step-2 recipe.
    Single linkage chains, complete linkage is the strictest; average sits
    between and is the usual default for expression profiles.
    """
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist
    X = np.nan_to_num(matrix.to_numpy(dtype=float).T, nan=0.0)
    labels = list(matrix.columns)
    if X.shape[0] < 2:
        return None, labels
    d = pdist(X, metric=metric)
    d = np.nan_to_num(d, nan=float(np.nanmax(d)) if np.isfinite(d).any() else 1.0)
    return linkage(d, method=method), labels


def cut_tree_clusters(Z, labels: list[str], ks: list[int]) -> pd.DataFrame:
    """Cluster membership at several cut heights, one column per k."""
    from scipy.cluster.hierarchy import fcluster
    out = pd.DataFrame({"label": labels})
    for k in ks:
        if 2 <= k <= len(labels):
            out[f"cluster_k{k}"] = fcluster(Z, t=k, criterion="maxclust")
    return out


def kmeans_clusters(matrix: pd.DataFrame, k: int, seed: int = 0, n_init: int = 10):
    """K-means over the ROWS (genes), the PDF's step-4 expression patterns.

    Lloyd's algorithm with k-means++ starts, best of `n_init` by inertia.
    Deterministic for a fixed seed so the exported clusters are reproducible.
    """
    from scipy.cluster.vq import kmeans2
    X = np.nan_to_num(matrix.to_numpy(dtype=float), nan=0.0)
    if X.shape[0] < k or k < 2:
        return np.zeros(X.shape[0], dtype=int), np.empty((0, X.shape[1])), np.nan
    best = (None, None, np.inf)
    for i in range(n_init):
        try:
            cent, lab = kmeans2(X, k, minit="++", seed=seed + i, missing="warn")
        except Exception:
            cent, lab = kmeans2(X, k, minit="points", seed=seed + i, missing="warn")
        if len(np.unique(lab)) < 2:
            continue
        inertia = float(((X - cent[lab]) ** 2).sum())
        if inertia < best[2]:
            best = (lab, cent, inertia)
    if best[0] is None:
        return np.zeros(X.shape[0], dtype=int), np.empty((0, X.shape[1])), np.nan
    return best[0], best[1], best[2]


def cluster_quality(matrix: pd.DataFrame, labels: np.ndarray) -> pd.DataFrame:
    """Per-cluster homogeneity and separation -- the trade-off the PDF stresses.

    Homogeneity is the mean correlation of a member to its cluster centroid
    (higher = tighter). Separation is the max correlation of this centroid to
    any other centroid (lower = more distinct). Raising k buys homogeneity and
    costs separation, and both numbers are needed to see that happening.
    """
    X = np.nan_to_num(matrix.to_numpy(dtype=float), nan=0.0)
    labels = np.asarray(labels)
    uniq = np.unique(labels)
    cents = np.vstack([X[labels == u].mean(axis=0) for u in uniq])

    def corr(a, b):
        a, b = a - a.mean(), b - b.mean()
        d = np.linalg.norm(a) * np.linalg.norm(b)
        return float(a @ b / d) if d > 0 else np.nan

    rows = []
    for i, u in enumerate(uniq):
        members = X[labels == u]
        homog = [corr(m, cents[i]) for m in members]
        others = [corr(cents[i], cents[j]) for j in range(len(uniq)) if j != i]
        rows.append({
            "cluster": int(u), "n_members": int(members.shape[0]),
            "homogeneity_mean_corr_to_centroid": float(np.nanmean(homog)) if homog else np.nan,
            "homogeneity_min": float(np.nanmin(homog)) if homog else np.nan,
            "separation_max_corr_to_other_centroid": float(np.nanmax(others)) if others else np.nan,
            "within_cluster_ss": float(((members - cents[i]) ** 2).sum()),
        })
    return pd.DataFrame(rows)


def silhouette_corr(matrix: pd.DataFrame, labels: np.ndarray) -> float:
    """Mean silhouette width using correlation distance (1 - Pearson r).

    For each gene: a = mean distance to its own cluster, b = mean distance to
    the nearest other cluster, s = (b - a) / max(a, b). Unlike homogeneity
    minus the max between-centroid correlation, this does not drift with K by
    construction, so it can pick K > 2 when the data hold more patterns.
    Singleton clusters score 0, the usual convention.
    """
    from scipy.spatial.distance import pdist, squareform
    X = np.nan_to_num(matrix.to_numpy(dtype=float), nan=0.0)
    labels = np.asarray(labels)
    uniq, idx = np.unique(labels, return_inverse=True)
    if len(uniq) < 2:
        return np.nan
    D = np.nan_to_num(squareform(pdist(X, metric="correlation")), nan=1.0)
    onehot = np.zeros((len(labels), len(uniq)))
    onehot[np.arange(len(labels)), idx] = 1.0
    counts = onehot.sum(axis=0)
    sums = D @ onehot                              # gene x cluster distance sums
    own_n = counts[idx]
    a = np.divide(sums[np.arange(len(labels)), idx], own_n - 1,
                  out=np.zeros(len(labels)), where=own_n > 1)
    mean_other = sums / counts
    mean_other[np.arange(len(labels)), idx] = np.inf
    b = mean_other.min(axis=1)
    s = np.where(own_n > 1, (b - a) / np.maximum(np.maximum(a, b), 1e-12), 0.0)
    return float(s.mean())


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    """Standardise each gene to mean 0 / SD 1 across groups.

    Required before clustering on *pattern*: without it a highly expressed gene
    and a lowly expressed one with identical kinetics land in different
    clusters purely because of their magnitude.
    """
    X = matrix.to_numpy(dtype=float)
    mu = np.nanmean(X, axis=1, keepdims=True)
    sd = np.nanstd(X, axis=1, keepdims=True)
    Z = np.divide(X - mu, sd, out=np.zeros_like(X), where=sd > 0)
    return pd.DataFrame(Z, index=matrix.index, columns=matrix.columns)


# ---------------------------------------------------------------------------
# GSEA
# ---------------------------------------------------------------------------
def gsea_enrichment_score(ranked_genes: list[str], gene_set: set[str],
                          weights: np.ndarray | None = None):
    """Weighted Kolmogorov-Smirnov running sum (Subramanian et al. 2005).

    Walks the ranked list adding weight for set members and subtracting for
    non-members; ES is the maximum deviation from zero. Unlike a hypergeometric
    test on a thresholded list, this uses the whole ranking, so it still detects
    a set that shifts coherently without any member clearing a cutoff -- the
    low-power problem the PDF raises for enrichment over DE clusters.

    Returns (ES, running_sum, leading_edge_genes).
    """
    n = len(ranked_genes)
    if n == 0:
        return np.nan, np.array([]), []
    in_set = np.array([g in gene_set for g in ranked_genes])
    n_hits = int(in_set.sum())
    if n_hits == 0 or n_hits == n:
        return np.nan, np.zeros(n), []

    w = np.ones(n) if weights is None else np.abs(np.asarray(weights, dtype=float))
    hit_w = np.where(in_set, w, 0.0)
    total = hit_w.sum()
    if total <= 0:
        hit_w, total = in_set.astype(float), float(n_hits)
    p_hit = np.cumsum(hit_w) / total
    p_miss = np.cumsum(~in_set) / (n - n_hits)
    running = p_hit - p_miss
    peak = int(np.argmax(np.abs(running)))
    es = float(running[peak])
    # Leading edge: set members up to (or beyond) the peak, by its sign.
    idx = np.arange(n)
    core = idx[: peak + 1] if es >= 0 else idx[peak:]
    leading = [ranked_genes[i] for i in core if in_set[i]]
    return es, running, leading


def _null_es(n: int, size: int, weights, n_permutations: int, seed: int,
             chunk: int = 250) -> np.ndarray:
    """ES of `n_permutations` random sets of `size` genes, vectorised.

    Draws the same random sets as calling gsea_enrichment_score once per
    permutation with rng.choice(universe, size, replace=False), so results are
    identical to the per-permutation loop -- just without the Python overhead.
    """
    rng = np.random.default_rng(seed)
    w = np.ones(n) if weights is None else np.abs(np.asarray(weights, dtype=float))
    out = np.empty(n_permutations)
    for start in range(0, n_permutations, chunk):
        m = min(chunk, n_permutations - start)
        hits = np.zeros((m, n), dtype=bool)
        for r in range(m):
            hits[r, rng.choice(n, size=size, replace=False)] = True
        hit_w = np.where(hits, w, 0.0)
        total = hit_w.sum(axis=1, keepdims=True)
        # Sets whose members all carry zero weight fall back to unweighted.
        zero = (total <= 0).ravel()
        if zero.any():
            hit_w[zero] = hits[zero]
            total[zero] = size
        running = np.cumsum(hit_w, axis=1) / total - np.cumsum(~hits, axis=1) / (n - size)
        peak = np.argmax(np.abs(running), axis=1)
        out[start:start + m] = running[np.arange(m), peak]
    return out


def gsea_test(ranked_genes: list[str], gene_set: set[str], weights=None,
              n_permutations: int = 1000, seed: int = 0,
              null_cache: dict | None = None) -> dict:
    """GSEA with a gene-set permutation null, as described in the course notes.

    Random sets of the same size give the null distribution of ES. Note the
    caveat: permuting gene *sets* (not sample labels) does not preserve
    gene-gene correlation, so the p-value is anti-conservative for correlated
    sets. It ranks hypotheses well; treat the absolute value with care.

    The null depends only on (ranking weights, set size, seed), so callers
    testing many sets against one ranking can pass a dict as `null_cache` and
    reuse it; the cache must be fresh for every new ranking.
    """
    es, running, leading = gsea_enrichment_score(ranked_genes, gene_set, weights)
    if not np.isfinite(es):
        return {"es": np.nan, "nes": np.nan, "p_value": np.nan,
                "n_genes_in_set": 0, "leading_edge": [], "peak_rank": np.nan}
    n = len(ranked_genes)
    size = sum(1 for g in ranked_genes if g in gene_set)
    if null_cache is not None and size in null_cache:
        null = null_cache[size]
    else:
        null = _null_es(n, size, weights, n_permutations, seed)
        null = null[np.isfinite(null)]
        if null_cache is not None:
            null_cache[size] = null
    if null.size == 0:
        return {"es": es, "nes": np.nan, "p_value": np.nan,
                "n_genes_in_set": size, "leading_edge": leading,
                "peak_rank": int(np.argmax(np.abs(running))) + 1}
    same = null[null >= 0] if es >= 0 else null[null < 0]
    # +1 numerator keeps p strictly positive (Phipson & Smyth).
    p = (np.sum(np.abs(same) >= abs(es)) + 1) / (same.size + 1) if same.size else np.nan
    nes = es / np.abs(np.mean(same)) if same.size and np.mean(same) != 0 else np.nan
    return {"es": es, "nes": float(nes) if np.isfinite(nes) else np.nan,
            "p_value": float(p), "n_genes_in_set": size,
            "leading_edge": leading,
            "peak_rank": int(np.argmax(np.abs(running))) + 1,
            "leading_edge_size": len(leading)}


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """BH-adjusted p-values (q-values) -- the PDF's FDR correction."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    n = pv.size
    q = pv[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
    adj = np.empty(n)
    adj[order] = q
    out[ok] = adj
    return out
