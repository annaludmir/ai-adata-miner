"""Shared helpers for step-3 analyses: CSV access, the replicate structure, stats.

Everything here reads csv_exports/ only -- never an h5ad -- so step 3 runs on a
laptop. Dependencies: numpy, pandas, scipy (matplotlib optional, figures only).

Design facts every analysis must respect (established by 01_data_audit):
  * Chemistry strata (<dataset>__v2 / __v3) have disjoint donors and
    near-disjoint ages, so they are analysed separately and used as two
    independent replicates of each other.
  * Within a stratum, almost every age comes from one donor: the replicate
    unit for anything about age is the donor (or the age point), never a cell.
  * cortex donors are a subset of human_dev donors. Agreement between the two
    datasets is reproducibility of processing, not independent replication.
"""
from __future__ import annotations

import itertools
import math
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
EXPORTS = Path(os.environ.get("AI_ADATA_OUT_ROOT", REPO / "csv_exports"))
RESULTS = Path(os.environ.get("AIM_DOWNSTREAM_OUT", REPO / "downstream_analyses" / "results"))
# Same rules file as extraction (config.EXCLUSIONS_FILE); 'none' disables.
EXCLUSIONS = os.environ.get("AIM_EXCLUSIONS", str(REPO / "exclusions.csv"))

DATASETS = ("cortex", "human_dev")
CHEMISTRIES = ("v2", "v3")
STRATA = [(d, c) for d in DATASETS for c in CHEMISTRIES]
SEED = 0

# Groups smaller than this are too noisy to anchor a mean expression profile.
MIN_CELLS = 50
# A donor contributing fewer cells than this is dropped from composition work.
MIN_DONOR_CELLS = 200

# Sex-linked genes used to infer donor sex (obs['sex'] is "unknown") and to
# flag findings that may be donor-sex rather than biology.
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY",
           "NLGN4Y", "TXLNGY", "PRKY"]
X_INACTIVATION_GENES = ["XIST", "TSIX"]

MARKER_PANEL_TO_CLASS = {
    "radial_glia": "Radial glia", "neuronal_ipc": "Neuronal IPC",
    "neuroblast": "Neuroblast", "neuron": "Neuron",
    "glioblast_opc": "Glioblast", "oligo": "Oligo",
    "immune_microglia": "Immune", "vascular_endothelial": "Vascular",
    "erythrocyte": "Erythrocyte", "neural_crest": "Neural crest",
}


# ---------------------------------------------------------------------------
# logging / IO
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def ns(dataset: str, chem: str | None = None) -> str:
    return dataset if chem is None else f"{dataset}__{chem}"


def require(path: Path) -> Path:
    """Fail loudly when an input is missing -- never analyse a guessed file."""
    if not path.exists():
        sys.exit(f"ERROR: required input missing: {path}\n"
                 f"  Run the extraction step that writes it (see the manifest), "
                 f"or check AI_ADATA_OUT_ROOT (now {EXPORTS}).")
    return path


def csv(namespace: str, rel: str, **kw) -> pd.DataFrame:
    return pd.read_csv(require(EXPORTS / namespace / rel), **kw)


def manifest(namespace: str) -> pd.DataFrame:
    return csv(namespace, "_manifest.csv")


class Output:
    """Result folder for one analysis: CSVs, figures and the SUMMARY.md."""

    def __init__(self, slug: str):
        check_numerics()
        self.slug = slug
        self.dir = RESULTS / slug
        self.dir.mkdir(parents=True, exist_ok=True)
        # Start clean: a figure or table the analysis no longer produces (e.g. a
        # scope that fell below the donor minimum) must not linger as a result.
        for old in self.dir.iterdir():
            if old.is_file():
                old.unlink()
        self.inputs: set[str] = set()
        self.written: list[tuple[str, str]] = []

    def used(self, *rels: str) -> None:
        self.inputs.update(rels)

    def write(self, df: pd.DataFrame, name: str, description: str,
              index: bool = False) -> Path:
        p = self.dir / f"{name}.csv"
        df.to_csv(p, index=index)
        self.written.append((p.name, description))
        log(f"  wrote {p.name} ({len(df):,} rows) -- {description}")
        return p

    def figure(self, fig, name: str, description: str) -> None:
        p = self.dir / f"{name}.png"
        fig.savefig(p, dpi=130, bbox_inches="tight")
        self.written.append((p.name, description))
        log(f"  wrote {p.name} -- {description}")

    def summary(self, title: str, question: str, method: list[str],
                findings: list[str], limitations: list[str],
                next_steps: list[str], extra: str = "") -> None:
        """SUMMARY.md in the shape the step-3 prompt asks for."""
        lines = [f"# {title}", "",
                 f"_Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} "
                 f"by `downstream_analyses/{self.slug}.py` from `{EXPORTS.name}/`._", "",
                 "## Question", "", question, "",
                 "## Inputs", ""]
        lines += [f"- `{r}`" for r in sorted(self.inputs)] or ["- (none)"]
        lines += ["", "## Method", ""] + [f"- {m}" for m in method]
        lines += ["", "## Key findings", ""] + [f"- {f}" for f in findings]
        if extra:
            lines += ["", extra.rstrip()]
        lines += ["", "## Limitations", ""] + [f"- {l}" for l in limitations]
        lines += ["", "## What would strengthen this", ""] + [f"- {n}" for n in next_steps]
        lines += ["", "## Output files", ""]
        lines += [f"- `{n}` -- {d}" for n, d in self.written]
        (self.dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")
        log(f"  wrote SUMMARY.md -> {self.dir}")


def plt_or_none():
    """matplotlib.pyplot, or None (figures are skipped, tables still written)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        log("  matplotlib not installed -- figures skipped")
        return None


# ---------------------------------------------------------------------------
# pseudobulk tables
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# exclusions -- re-applied here so exports made before a rule was added are
# still filtered wherever a table resolves age, donor or sample
# ---------------------------------------------------------------------------
_EXCL_CACHE: dict[str, pd.DataFrame] = {}


def exclusion_rules(dataset: str) -> pd.DataFrame:
    if dataset not in _EXCL_CACHE:
        cols = ["dataset", "role", "value", "reason"]
        if EXCLUSIONS.lower() == "none" or not Path(EXCLUSIONS).exists():
            rules = pd.DataFrame(columns=cols)
        else:
            rules = pd.read_csv(EXCLUSIONS, dtype=str, comment="#").fillna("")
        rules = rules[rules["dataset"].isin([dataset, "*"])]
        other = rules[~rules["role"].isin(["age", "donor", "sample"])]
        if len(other):
            log(f"  note: exclusion roles {sorted(set(other.role))} apply at extraction "
                "only; step 3 can re-apply age, donor and sample rules")
        _EXCL_CACHE[dataset] = rules.reset_index(drop=True)
    return _EXCL_CACHE[dataset]


def excluded(dataset: str, age=None, donor=None, sample=None) -> bool:
    """True if any exclusion rule for this dataset matches the given labels."""
    for r in exclusion_rules(dataset).itertuples():
        if r.role == "age" and age is not None and abs(float(age) - float(r.value)) < 0.01:
            return True
        if r.role == "donor" and donor is not None and normalise_donor(donor) == normalise_donor(r.value):
            return True
        if r.role == "sample" and sample is not None and str(sample) == r.value:
            return True
    return False


def _group_excluded(namespace: str, grouping: str) -> set[str]:
    ds = namespace.split("__")[0]
    if exclusion_rules(ds).empty:
        return set()
    gs = group_summary(namespace, grouping)
    labels = gs["group"].astype(str)
    if grouping == "age":
        return {g for g in labels if excluded(ds, age=float(g))}
    if grouping.endswith("_x_age"):
        return {g for g in labels if excluded(ds, age=split_class_age(g)[1])}
    if grouping == "donor":
        ages = _raw_donor_ages(namespace)
        return {g for g in labels if excluded(ds, donor=g, age=ages.get(g))}
    if grouping == "sample":
        meta = csv(namespace, "01_overview/sample_summary.csv").set_index("sample")
        return {g for g in labels if g in meta.index and excluded(
            ds, sample=g, donor=meta.loc[g, "donor"], age=meta.loc[g, "age_pcw"])}
    return set()


def group_summary(namespace: str, grouping: str) -> pd.DataFrame:
    return csv(namespace, f"09_pseudobulk/{grouping}__group_summary.csv")


def group_matrix(namespace: str, grouping: str, kind: str = "mean_lognorm",
                 min_cells: int = MIN_CELLS) -> pd.DataFrame:
    """gene x group table from 09_pseudobulk, keeping groups with >= min_cells.

    The n_cells filter also removes the all-zero columns of groups that exist
    only in the other chemistry (present in exports made before that fix).
    """
    df = csv(namespace, f"09_pseudobulk/{grouping}__{kind}.csv")
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    df.index.name = "gene"
    gs = group_summary(namespace, grouping)
    good = set(gs.loc[gs["n_cells"] >= min_cells, "group"].astype(str))
    good -= _group_excluded(namespace, grouping)
    return df[[c for c in df.columns if str(c) in good]].apply(pd.to_numeric, errors="coerce")


def tmm_log_cpm(counts: pd.DataFrame, prior: float = 1.0) -> pd.DataFrame:
    """log2(TMM-normalised CPM + prior) from a gene x group raw-count table.

    Use this, not mean_lognorm, when comparing groups that differ in depth or in
    transcriptome composition. mean_lognorm averages log1p(CP10K) over cells, so
    it drifts with UMIs per cell; and a few dominant genes (haemoglobin in
    erythrocytes) depress every other gene's CPM. TMM corrects both from the
    bulk of unchanged genes (Robinson & Oshlack 2010; the repo's validated
    implementation in lib/bulk_stats.py).
    """
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from lib.bulk_stats import tmm_normalization_factors
    tmm = tmm_normalization_factors(counts).set_index("sample")
    eff = (tmm["library_size"] * tmm["tmm_factor"]).reindex(counts.columns)
    cpm = counts.to_numpy(float) / eff.to_numpy()[None, :] * 1e6
    return pd.DataFrame(np.log2(cpm + prior), index=counts.index, columns=counts.columns)


def split_class_age(label: str) -> tuple[str, float]:
    """'Radial glia | 6.9' -> ('Radial glia', 6.9)."""
    cls, age = label.rsplit("|", 1)
    return cls.strip(), float(age)


# ---------------------------------------------------------------------------
# donors, ages, panels, genes
# ---------------------------------------------------------------------------
def normalise_donor(d: str) -> str:
    """cortex writes 'XHU:1966:307' where human_dev writes 'XHU:307'."""
    return re.sub(r"^(XHU|XDD):\d+:(\d+)$", r"\1:\2", str(d))


def _raw_donor_ages(namespace: str) -> dict[str, float]:
    x = csv(namespace, "05_confounds/crosstab_age_x_donor.csv", index_col=0)
    return {d: float(x[d].idxmax()) for d in x.columns if x[d].sum() > 0}


def donor_ages(namespace: str) -> pd.DataFrame:
    """donor, age_pcw, n_cells for one stratum (each donor has exactly one age).

    Donors matched by an exclusion rule (by donor or by age) are left out.
    """
    ds = namespace.split("__")[0]
    x = csv(namespace, "05_confounds/crosstab_age_x_donor.csv", index_col=0)
    rows = []
    for donor in x.columns:
        col = x[donor]
        if col.sum() > 0 and not excluded(ds, donor=donor, age=float(col.idxmax())):
            rows.append({"donor": donor, "age_pcw": float(col.idxmax()),
                         "n_ages": int((col > 0).sum()), "n_cells": int(col.sum())})
    return pd.DataFrame(rows)


def panels(namespace: str) -> pd.DataFrame:
    """panel_group, panel, gene, present_in_dataset, exported_in_pseudobulk."""
    return csv(namespace, "11_panels/panel_coverage.csv")


def chromosome_map() -> dict[str, str]:
    """gene symbol -> chromosome ('X', 'Y', '1', ...), from the cross-dataset map."""
    g = csv("_cross_dataset", "gene_id_map.csv", low_memory=False)
    g = g.dropna(subset=["symbol", "Chromosome"])
    chrom = g["Chromosome"].astype(str).str.replace("^chr", "", regex=True)
    return dict(zip(g["symbol"].astype(str), chrom))


def sex_linked(genes, chrom: dict[str, str]) -> np.ndarray:
    """True for Y-chromosome genes and XIST/TSIX: donor-sex signal, not biology."""
    return np.array([chrom.get(g) == "Y" or g in X_INACTIVATION_GENES for g in genes])


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------
def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * pv.size / (np.arange(pv.size) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
        out = np.empty(pv.size)
        out[order] = ranked
        q[ok] = out
    return q


# ---------------------------------------------------------------------------
# numerical gate: refuse to run on a numpy build that computes wrongly
# ---------------------------------------------------------------------------
# numpy 2.2.6 on macOS 26.2 (the Mac's system Python) silently computed wrong
# rank correlations on large arrays -- 8,000 x 6 matrices gave Spearman rho
# above 1 -- with the error moving between operations (ranking, centring,
# products) depending on how arrays sat in memory, so checking any single
# operation was not enough. numpy 2.5.3 on the same machine and the cluster's
# Linux numpy were correct. Every analysis therefore starts by running the
# real correlation routine on large inputs with known answers and stops with
# instructions if anything is off.
_NUMERICS_OK = False


def _numerics_fail(what: str) -> None:
    sys.exit("ERROR: numpy " + np.__version__ + " on this machine computes wrong results (" + what
             + "). Seen with numpy 2.2.6 on macOS 26; every statistic would be silently wrong. "
             "Use a newer numpy in a separate environment, e.g.\n"
             "  python3 -m venv ~/aim-env && ~/aim-env/bin/pip install -U numpy scipy pandas matplotlib\n"
             "  PYTHON=~/aim-env/bin/python ./downstream_analyses/run_all.sh\n"
             "or run step 3 on the cluster (slurm_04_downstream.sh).")


def check_numerics() -> None:
    """End-to-end check of the correlation machinery on large inputs with known answers."""
    global _NUMERICS_OK
    if _NUMERICS_OK:
        return
    ages = np.array([7.5, 8.0, 8.5, 9.2, 9.5, 10.0])
    base = np.array([[-2.899, -3.711, -3.135, -2.303, -3.062, -3.332],
                     [5.786, 4.733, 4.854, 5.623, 5.572, 5.673]])
    known = np.array([stats.spearmanr(r, ages).statistic for r in base])
    try:
        for reps in (2000, 4000, 8000, 16000):
            rho = spearman_rows(np.vstack([base] * reps), ages)
            if not np.allclose(rho, np.tile(known, reps)):
                raise ValueError(f"rank correlation on {2 * reps} x 6 rows")
        rng = np.random.default_rng(7)
        X = rng.normal(size=(20000, 7)).round(3)
        x = np.arange(7.0)
        rho = spearman_rows(X, x)
        for i in rng.choice(len(X), size=200, replace=False):
            if not np.isclose(rho[i], stats.spearmanr(X[i], x).statistic):
                raise ValueError("rank correlation on 20,000 x 7 random rows")
        A, B = rng.normal(size=(3000, 400)), rng.normal(size=(400, 300))
        P = A @ B
        for _ in range(100):
            i, j = int(rng.integers(3000)), int(rng.integers(300))
            if not np.isclose(P[i, j], math.fsum(float(u) * float(v) for u, v in zip(A[i], B[:, j]))):
                raise ValueError("matrix product 3000 x 400 x 300")
    except ValueError as err:
        _numerics_fail(str(err))
    _NUMERICS_OK = True


def dot(a, b) -> np.ndarray:
    """a @ b for dense arrays (the single place products go)."""
    return np.asarray(a, dtype=float) @ np.asarray(b, dtype=float)


# Second line of defence, on every call: a correlation outside [-1, 1], or a
# few sampled entries that disagree with a pure-Python recomputation, stop the
# run. Cheap (a handful of short loops per call) and independent of numpy.
_SPOT_CHECKS = 6
_spot_rng = np.random.default_rng(SEED)


def _pearson_py(a, b) -> float:
    a, b = [float(u) for u in a], [float(u) for u in b]
    ma, mb = math.fsum(a) / len(a), math.fsum(b) / len(b)
    da, db = [u - ma for u in a], [u - mb for u in b]
    den = math.sqrt(math.fsum(u * u for u in da) * math.fsum(u * u for u in db))
    return math.fsum(u * w for u, w in zip(da, db)) / den if den > 0 else 0.0


def _rank_py(v) -> list[float]:
    """Average ranks (ties share the mean rank), as scipy.stats.rankdata."""
    v = [float(u) for u in v]
    order = sorted(range(len(v)), key=v.__getitem__)
    ranks = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def _spot_check(values: np.ndarray, reference, what: str) -> None:
    """values: correlations just computed; reference(index) recomputes one in pure Python."""
    if np.nanmax(np.abs(values), initial=0.0) > 1 + 1e-9:
        _numerics_fail(f"a {what} above 1 in absolute value")
    ok = np.flatnonzero(np.isfinite(values.ravel()))
    for flat in _spot_rng.choice(ok, size=min(_SPOT_CHECKS, ok.size), replace=False):
        idx = np.unravel_index(flat, values.shape)
        if not math.isclose(values[idx], reference(idx), rel_tol=1e-7, abs_tol=1e-9):
            _numerics_fail(f"{what} at {tuple(int(i) for i in idx)} of {values.shape}")


def corr_rows(A: np.ndarray, B: np.ndarray | None = None) -> np.ndarray:
    """Pearson correlation between the rows of A (and of B), via dot()."""
    def z(M):
        M = np.asarray(M, dtype=float)
        M = M - M.mean(axis=1, keepdims=True)
        n = np.linalg.norm(M, axis=1, keepdims=True)
        return np.divide(M, n, out=np.zeros_like(M), where=n > 0)
    A0 = np.asarray(A, dtype=float)
    B0 = A0 if B is None else np.asarray(B, dtype=float)
    A = z(A0)
    B = A if B is None else z(B0)
    out = dot(A, B.T)
    _spot_check(out, lambda ij: _pearson_py(A0[ij[0]], B0[ij[1]]), "Pearson correlation")
    return out


def rank_rows(X: np.ndarray) -> np.ndarray:
    return np.apply_along_axis(stats.rankdata, 1, X)


def spearman_rows(X: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Spearman rho of every row of X against the vector x (ties averaged)."""
    X = np.asarray(X, dtype=float)
    R = rank_rows(X)
    rx = stats.rankdata(x)
    R = R - R.mean(axis=1, keepdims=True)
    rx = rx - rx.mean()
    den = np.sqrt((R ** 2).sum(axis=1) * (rx ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.where(den > 0, dot(R, rx) / den, np.nan)
    _spot_check(rho, lambda i: _pearson_py(_rank_py(X[i[0]]), _rank_py(x)), "rank correlation")
    return rho


def permutation_orders(n: int, n_perm: int = 200_000, seed: int = SEED) -> np.ndarray:
    """All n! orderings when that is <= n_perm (exact test), else a random sample."""
    if math.factorial(n) <= n_perm:
        return np.array(list(itertools.permutations(range(n))), dtype=np.int32)
    rng = np.random.default_rng(seed)
    return np.argsort(rng.random((n_perm, n)), axis=1).astype(np.int32)


def spearman_perm_p(rho: np.ndarray, x: np.ndarray, X: np.ndarray | None = None,
                    n_perm: int = 200_000, seed: int = SEED) -> tuple[np.ndarray, bool]:
    """Two-sided permutation p for Spearman rho against x.

    With few points (5-11 ages or donors) the t approximation behind
    scipy.stats.spearmanr is unreliable, so the null is built by permuting x:
    exactly over all n! orders when feasible. When X (the rows rho came from)
    is None, rows are assumed tie-free and one shared null is used; otherwise
    each row gets its own null (correct with ties, slower -- use for few rows).
    Returns (p, exact).
    """
    x = np.asarray(x, dtype=float)
    orders = permutation_orders(len(x), n_perm, seed)
    exact = math.factorial(len(x)) <= n_perm
    rho = np.asarray(rho, dtype=float)
    if X is None:
        rx = stats.rankdata(x)
        base = np.arange(1, len(x) + 1, dtype=float)
        base -= base.mean()
        rxp = rx[orders] - rx.mean()
        null = dot(rxp, base) / np.sqrt((base ** 2).sum() * (rxp ** 2).sum(axis=1))
        null = np.sort(np.abs(null))
        # count |null| >= |rho| with a small tolerance for float ties
        k = null.size - np.searchsorted(null, np.abs(rho) - 1e-12, side="left")
        p = (k + (0 if exact else 1)) / (null.size + (0 if exact else 1))
        return np.where(np.isfinite(rho), p, np.nan), exact
    R = rank_rows(np.asarray(X, dtype=float))
    R = R - R.mean(axis=1, keepdims=True)
    rx = stats.rankdata(x)
    rxp = rx[orders] - rx.mean()
    den = np.sqrt((R ** 2).sum(axis=1)[:, None] * (rxp ** 2).sum(axis=1)[None, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        null = np.abs(np.where(den > 0, dot(R, rxp.T) / den, 0.0))
    k = (null >= np.abs(rho)[:, None] - 1e-12).sum(axis=1)
    p = (k + (0 if exact else 1)) / (orders.shape[0] + (0 if exact else 1))
    return np.where(np.isfinite(rho), p, np.nan), exact


def signed_stouffer(rhos: list[np.ndarray], ps: list[np.ndarray],
                    weights: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Combine two-sided tests from independent replicates, keeping direction.

    z_i = sign(rho_i) * Phi^-1(1 - p_i/2); Z = sum(w_i z_i) / sqrt(sum w_i^2).
    Agreement in direction strengthens Z; disagreement cancels it -- which is
    the point: a trend that flips between donor sets is not replicated.
    Returns (Z, two-sided p).
    """
    zs = []
    for r, p in zip(rhos, ps):
        p = np.clip(np.asarray(p, dtype=float), 1e-300, 1.0)
        zs.append(np.sign(r) * stats.norm.isf(p / 2))
    w = np.asarray(weights, dtype=float)
    Z = sum(wi * zi for wi, zi in zip(w, zs)) / np.sqrt((w ** 2).sum())
    return Z, 2 * stats.norm.sf(np.abs(Z))


def replication_tier(rho_a, p_a, rho_b, p_b, q, alpha: float = 0.05) -> np.ndarray:
    """'replicated' / 'supported' / '' for a trend tested in two donor sets.

    replicated: same direction, combined q < alpha, and each set nominally
                significant on its own (one-sided p < alpha, i.e. two-sided
                p < 2*alpha, in the shared direction).
    supported:  same direction and combined q < alpha, but one set is weak --
                the combined test is carried by the other.
    """
    rho_a, rho_b = np.asarray(rho_a, float), np.asarray(rho_b, float)
    p_a, p_b, q = (np.asarray(x, float) for x in (p_a, p_b, q))
    same = np.sign(rho_a) == np.sign(rho_b)
    comb = same & (q < alpha)
    both = (p_a < 2 * alpha) & (p_b < 2 * alpha)
    return np.where(comb & both, "replicated", np.where(comb, "supported", ""))


def clr(counts: pd.DataFrame, pseudo: float = 0.5) -> pd.DataFrame:
    """Centred log-ratio per row: composition-aware, unlike raw fractions."""
    lg = np.log(counts.to_numpy(dtype=float) + pseudo)
    return pd.DataFrame(lg - lg.mean(axis=1, keepdims=True),
                        index=counts.index, columns=counts.columns)


# ---------------------------------------------------------------------------
# per-(cell class, age) tables and cell-class preference of gene sets
# ---------------------------------------------------------------------------
def zscore_rows(X: np.ndarray) -> np.ndarray:
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True)
    return np.divide(X - mu, sd, out=np.zeros_like(X), where=sd > 0)


def class_age_logcpm(n: str):
    """{cell_class: (log2 TMM-CPM genes x ages, detection genes x ages, ages, min_cells)}.

    TMM is computed within a class, across its age points: the question is how
    one cell type changes, so its own ages are the libraries to align.
    """
    cnt = group_matrix(n, "cell_class_x_age", "pseudobulk_counts")
    det = group_matrix(n, "cell_class_x_age", "detection_fraction")
    gs = group_summary(n, "cell_class_x_age").set_index("group")
    out = {}
    by_class: dict[str, list[tuple[float, str]]] = {}
    for col in cnt.columns:
        cls, age = split_class_age(col)
        by_class.setdefault(cls, []).append((age, col))
    for cls, items in by_class.items():
        items.sort()
        cols = [c for _, c in items]
        ages = np.array([a for a, _ in items])
        out[cls] = (tmm_log_cpm(cnt[cols]), det[cols], ages,
                    int(gs.loc[cols, "n_cells"].min()))
    return out


def set_mean_rows(Z: np.ndarray, sets: np.ndarray) -> np.ndarray:
    """Mean of Z's rows over each set: (n_sets x k indices) -> (n_sets x columns).

    Same as Z[sets].mean(axis=1) but through a sparse sum, so 5,000 random sets
    of a 1,000-gene list do not materialise a sets x genes x columns array.
    Repeated indices count once per occurrence, as in the dense version.
    """
    from scipy import sparse
    sets = np.asarray(sets)
    n_sets, k = sets.shape
    ind = sparse.csr_matrix((np.ones(sets.size), (np.repeat(np.arange(n_sets), k), sets.ravel())),
                            shape=(n_sets, Z.shape[0]))
    return np.asarray(ind @ Z) / k


def set_class_preference(n: str, ds: str, chem: str, sets: dict[str, list[str]],
                         rng: np.random.Generator, n_random: int = 5000,
                         n_bins: int = 10, min_genes: int = 5, min_age_points: int = 3,
                         label: str = "panel", match_extra: pd.Series | None = None) -> pd.DataFrame:
    """Does each gene set score higher in one cell class than in the others?

    Pseudobulk counts per (cell class, age point) -> log2 TMM-CPM; genes >= 5
    CPM on average; each gene Z-scored across all columns; a set's score in a
    column is the mean Z of its genes. For class c, at every age point, the
    difference between c and the mean of the other classes is taken, and T is
    its mean over age points (age points are donors). The null is T for random
    sets drawing one gene from each member's expression decile (and, with
    match_extra -- a per-gene category indexed by symbol, e.g. a length
    tertile -- from the same category too).
    """
    cnt = group_matrix(n, "cell_class_x_age", "pseudobulk_counts")
    lc = tmm_log_cpm(cnt)
    lc = lc.loc[lc.mean(axis=1) >= np.log2(5 + 1)]
    labels = [split_class_age(c) for c in lc.columns]
    cls_of = np.array([c for c, _ in labels])
    age_of = np.array([a for _, a in labels])
    Z = zscore_rows(lc.to_numpy(float))
    gidx = {g: i for i, g in enumerate(lc.index.to_numpy())}

    level = lc.mean(axis=1).to_numpy()
    bins = np.digitize(level, np.quantile(level, np.linspace(0, 1, n_bins + 1)[1:-1]))
    if match_extra is not None:
        codes = pd.factorize(pd.Series(match_extra).reindex(lc.index), use_na_sentinel=True)[0]
        bins = bins * (codes.max() + 2) + (codes + 1)
    members_of_bin = {b: np.nonzero(bins == b)[0] for b in np.unique(bins)}

    classes = sorted(set(cls_of))
    layout = {}
    for c in classes:
        pairs = []
        for a in sorted(set(age_of)):
            here = np.nonzero(age_of == a)[0]
            mine = [j for j in here if cls_of[j] == c]
            others = [j for j in here if cls_of[j] != c]
            if mine and others:
                pairs.append((mine[0], others))
        layout[c] = pairs

    def T(score_cols: np.ndarray, c: str) -> tuple[np.ndarray, np.ndarray]:
        d = np.stack([score_cols[:, m] - score_cols[:, o].mean(axis=1)
                      for m, o in layout[c]], axis=1)
        return d.mean(axis=1), (d > 0).mean(axis=1)

    rows = []
    for name, genes in sets.items():
        idx = np.array([gidx[x] for x in genes if x in gidx])
        if idx.size < min_genes:
            continue
        obs = Z[idx].mean(axis=0)[None, :]
        rand = np.empty((n_random, idx.size), dtype=int)
        for k, gi in enumerate(idx):
            rand[:, k] = rng.choice(members_of_bin[bins[gi]], size=n_random)
        null_scores = set_mean_rows(Z, rand)
        for c in classes:
            if len(layout[c]) < min_age_points:
                continue
            t_obs, frac = T(obs, c)
            t_null, _ = T(null_scores, c)
            t_obs, frac = float(t_obs[0]), float(frac[0])
            p = (np.sum(np.abs(t_null - t_null.mean()) >= abs(t_obs - t_null.mean())) + 1) / (n_random + 1)
            rows.append({"dataset": ds, "chemistry": chem, label: name, "cell_class": c,
                         "n_genes": int(idx.size), "n_age_points": len(layout[c]),
                         "T_mean_z_difference": t_obs, "null_mean": float(t_null.mean()),
                         "null_sd": float(t_null.std()),
                         "effect_vs_null_sd": (t_obs - t_null.mean()) / t_null.std()
                         if t_null.std() > 0 else np.nan,
                         "frac_age_points_higher": frac, "perm_p": p})
    return pd.DataFrame(rows)


def combine_preference(pref: pd.DataFrame, label: str = "panel") -> pd.DataFrame:
    """v2 x v3 per (dataset, set, class): signed Stouffer, BH per dataset, tier."""
    comb = []
    for (ds, name, c), g in pref.groupby(["dataset", label, "cell_class"], sort=False):
        g = g.set_index("chemistry")
        if not set(CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        eff = [a.effect_vs_null_sd, b.effect_vs_null_sd]
        Z, pc = signed_stouffer([np.array([eff[0]]), np.array([eff[1]])],
                                [np.array([a.perm_p]), np.array([b.perm_p])],
                                [np.sqrt(a.n_age_points), np.sqrt(b.n_age_points)])
        comb.append({"dataset": ds, label: name, "cell_class": c,
                     "effect_v2": eff[0], "p_v2": a.perm_p, "frac_higher_v2": a.frac_age_points_higher,
                     "effect_v3": eff[1], "p_v3": b.perm_p, "frac_higher_v3": b.frac_age_points_higher,
                     "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
    comb = pd.DataFrame(comb)
    if comb.empty:
        return comb
    comb["combined_q"] = np.nan
    for _, idx in comb.groupby("dataset").groups.items():
        comb.loc[idx, "combined_q"] = bh(comb.loc[idx, "combined_p"])
    comb["tier"] = replication_tier(comb.effect_v2, comb.p_v2, comb.effect_v3, comb.p_v3,
                                    comb.combined_q)
    comb["direction"] = np.where(comb.stouffer_z > 0, "enriched", "depleted")
    return comb.sort_values(["dataset", "combined_p"])


# ---------------------------------------------------------------------------
# user gene lists
# ---------------------------------------------------------------------------
def gene_lists_dir() -> Path:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    import config
    return Path(config.GENE_LISTS_DIR)


_LISTS_CACHE: dict[str, dict[str, list[str]]] = {}
_MAPPED_CACHE: dict[tuple, dict[str, list[str]]] = {}


def gene_lists_label() -> str:
    """The gene-list folder for reports: repo-relative when inside the repo."""
    d = gene_lists_dir()
    try:
        return str(d.resolve().relative_to(REPO.resolve()))
    except ValueError:
        return str(d)


def gene_lists() -> dict[str, list[str]]:
    """{list name: genes as written in the file}, from config.GENE_LISTS_DIR.

    Same reader as the extraction step (lib/panels.read_gene_list): a CSV with a
    'gene' column (or its first column), or plain text with one gene per line.
    Set AIM_GENE_LISTS to point elsewhere.
    """
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from lib.panels import load_user_lists
    key = str(gene_lists_dir())
    if key not in _LISTS_CACHE:
        _LISTS_CACHE[key] = load_user_lists(gene_lists_dir())
    return _LISTS_CACHE[key]


_ID_MAP: dict[str, pd.DataFrame] = {}


def map_genes(genes, dataset: str) -> pd.DataFrame:
    """Resolve list entries to this dataset's gene symbols.

    In order: exact symbol; case-insensitive; Ensembl id (version stripped);
    the OTHER file's annotation -- a symbol known there is carried over through
    the shared Ensembl id; and, when the HGNC table has been fetched
    (fetch_annotations.sh), HGNC previous symbols and unambiguous aliases,
    through their Ensembl id. The two files use different symbol
    versions (cortex has HIST1H1C where human_dev has H1-2 after the 2020 HGNC
    histone renaming), so this recovers renamed genes without an alias table.
    Returns input, symbol and how it matched ('missing' if not).
    """
    if "_all" not in _ID_MAP:
        _ID_MAP["_all"] = csv("_cross_dataset", "gene_id_map.csv", low_memory=False)
    allg = _ID_MAP["_all"]
    g = allg[allg["dataset"] == dataset]
    o = allg[allg["dataset"] != dataset]
    symbols = set(g["symbol"].astype(str))
    upper = {}
    for x in g["symbol"].astype(str):
        upper.setdefault(x.upper(), x)
    ensg = dict(zip(g["accession_base"].astype(str), g["symbol"].astype(str)))
    other = {}
    for sym, acc in zip(o["symbol"].astype(str), o["accession_base"].astype(str)):
        other.setdefault(sym.upper(), acc)
    hgnc = hgnc_to_ensembl()     # empty until the HGNC table is fetched (part C)
    rows = []
    for x in genes:
        x = str(x).strip()
        if x in symbols:
            rows.append((x, x, "exact"))
        elif x.upper() in upper:
            rows.append((x, upper[x.upper()], "case"))
        elif x.split(".")[0] in ensg:
            rows.append((x, ensg[x.split(".")[0]], "ensembl"))
        elif x.upper() in other and other[x.upper()] in ensg:
            rows.append((x, ensg[other[x.upper()]], "other_file_symbol"))
        elif hgnc.get(x.upper(), "") in ensg:
            rows.append((x, ensg[hgnc[x.upper()]], "hgnc_previous_or_alias"))
        else:
            rows.append((x, "", "missing"))
    return pd.DataFrame(rows, columns=["input", "symbol", "match"])


# GWAS-derived lists name every gene near an associated variant, so one locus
# can contribute many co-regulated neighbours (16 HIST1 histones at 6p22 in a
# bipolar list made it look "mitotic"). Lists whose name matches this regex are
# collapsed to one gene per locus, represented by its most expressed member;
# 'none' disables. Curated lists are not collapsed.
COLLAPSE_PATTERN = os.environ.get("AIM_COLLAPSE_LISTS", "GWAS")
LOCUS_WINDOW = 1_000_000    # genes within this distance (bp) chain into one locus
_COORDS: dict[str, pd.DataFrame] = {}


def collapses(name: str) -> bool:
    return (COLLAPSE_PATTERN.lower() != "none"
            and re.search(COLLAPSE_PATTERN, name, flags=re.IGNORECASE) is not None)


_LEVEL: dict[str, pd.Series] = {}


def gene_level(dataset: str) -> pd.Series:
    """Mean CPM per gene over a dataset's cell classes and both chemistries."""
    if dataset not in _LEVEL:
        parts = []
        for chem in CHEMISTRIES:
            cnt = group_matrix(ns(dataset, chem), "cell_class", "pseudobulk_counts")
            parts.append((cnt.div(cnt.sum(axis=0), axis=1) * 1e6).mean(axis=1))
        _LEVEL[dataset] = pd.concat(parts, axis=1).mean(axis=1)
    return _LEVEL[dataset]


def gene_loci(symbols: list[str], dataset: str) -> pd.DataFrame:
    """Group genes into loci by position: same chromosome, chained within LOCUS_WINDOW.

    Returns symbol, chromosome, start, end, locus id, genes in the locus, mean
    CPM, and whether the gene represents its locus: the most highly expressed
    member -- the gene these data measure best. Position in the file is no
    guide: GWAS lists are typically sorted by coordinate, so "first listed"
    means "leftmost", often a non-coding gene the data cannot see. Every null
    in 06-08 matches on expression level, so this choice does not bias the
    tests. Genes without coordinates are their own locus.
    """
    if dataset not in _COORDS:
        g = csv("_cross_dataset", "gene_id_map.csv", low_memory=False)
        g = g[g["dataset"] == dataset].drop_duplicates("symbol")
        _COORDS[dataset] = pd.DataFrame({
            "chrom": g["Chromosome"].astype(str).str.replace("^chr", "", regex=True).to_numpy(),
            "start": pd.to_numeric(g["Start"], errors="coerce").to_numpy(),
            "end": pd.to_numeric(g["End"], errors="coerce").to_numpy()},
            index=g["symbol"].astype(str).to_numpy())
    co = _COORDS[dataset]
    df = pd.DataFrame({"symbol": symbols, "order": range(len(symbols))})
    df = df.join(co, on="symbol")
    df["locus"] = -1
    has = df["start"].notna() & df["end"].notna()
    next_id = 0
    for _, grp in df[has].sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        reach = -np.inf
        for i, r in grp.iterrows():
            if r.start > reach + LOCUS_WINDOW:
                next_id += 1
            df.loc[i, "locus"] = next_id
            reach = max(reach, r.end) if r.start <= reach + LOCUS_WINDOW else r.end
    for i in df.index[~has]:
        next_id += 1
        df.loc[i, "locus"] = next_id
    df["genes_in_locus"] = df.groupby("locus")["symbol"].transform("size")
    df["mean_cpm"] = df["symbol"].map(gene_level(dataset)).fillna(0.0)
    best = df.sort_values(["mean_cpm", "order"], ascending=[False, True]).drop_duplicates("locus").index
    df["representative"] = df.index.isin(best)
    return df.sort_values("order").drop(columns="order")


def mapped_lists(dataset: str, collapse: bool = True) -> dict[str, list[str]]:
    """User lists as unique symbols of this dataset (unmatched entries dropped).

    With collapse=True, lists matching COLLAPSE_PATTERN keep one gene per locus.
    """
    key = (dataset, collapse)
    if key not in _MAPPED_CACHE:
        out = {}
        for name, genes in gene_lists().items():
            m = map_genes(genes, dataset)
            syms = list(dict.fromkeys(m.loc[m.match != "missing", "symbol"]))
            if collapse and collapses(name):
                loci = gene_loci(syms, dataset)
                syms = list(loci.loc[loci.representative, "symbol"])
            out[name] = syms
        _MAPPED_CACHE[key] = out
    return _MAPPED_CACHE[key]


def analysis_gene_sets(dataset: str, groups=("ndd", "cell_cycle"),
                       modules: bool = True, out=None) -> dict[str, list[str]]:
    """Gene sets most analyses profile: user lists (GWAS lists collapsed per
    locus), seed panels of the given groups, and robust 08 modules."""
    sets = {f"list:{k}": v for k, v in mapped_lists(dataset).items()} if gene_lists() else {}
    n = ns(dataset, "v2")
    if out is not None:
        out.used(f"{n}/11_panels/panel_coverage.csv")
    pan = panels(n)
    for (_, pname), g in pan[pan.panel_group.isin(list(groups))].groupby(["panel_group", "panel"]):
        sets[f"seed:{pname}"] = list(g.gene)
    mod_path = RESULTS / "08_coexpression_modules" / "modules.csv"
    if modules and mod_path.exists():
        if out is not None:
            out.used("results/08_coexpression_modules/modules.csv")
        mods = pd.read_csv(mod_path)
        for r in mods[(mods.dataset == dataset) & mods.robust].itertuples():
            sets[f"module:{r.module}"] = r.genes.split("|")
    return sets


# ---------------------------------------------------------------------------
# co-expression across fine clusters
# ---------------------------------------------------------------------------
# The finest clustering in each file, and the minimum cells for a cluster's
# pseudobulk to be stable enough for a correlation.
COEXPR_GROUPING = {"cortex": "cluster_ClustersSurprise", "human_dev": "cluster_cluster_id"}
COEXPR_MIN_CELLS = 100
COEXPR_MIN_CPM = 10.0       # a gene must reach this CPM ...
COEXPR_MIN_CLUSTERS = 3     # ... in at least this many clusters


def cluster_expression(namespace: str):
    """(log2 TMM-CPM genes x clusters, cluster annotation) for co-expression.

    Clusters need >= COEXPR_MIN_CELLS cells and must not be donor-dominated
    (04's batch_dominated_flag: > 90% of cells from one donor), since a
    one-donor cluster can make a donor's quirks look like co-expression.
    Genes must reach COEXPR_MIN_CPM in >= COEXPR_MIN_CLUSTERS clusters.
    """
    ds = namespace.split("__")[0]
    grouping = COEXPR_GROUPING[ds]
    cnt = group_matrix(namespace, grouping, "pseudobulk_counts", min_cells=COEXPR_MIN_CELLS)
    clustering = grouping.replace("cluster_", "", 1)
    prof = csv(namespace, f"04_clusters/cluster_profile_{clustering}.csv")
    prof["cluster"] = prof["cluster"].astype(str)
    prof = prof.set_index("cluster")
    keep = [c for c in cnt.columns
            if c in prof.index and not bool(prof.loc[c, "batch_dominated_flag"])]
    cnt = cnt[keep]
    lc = tmm_log_cpm(cnt)
    expressed = (lc >= np.log2(COEXPR_MIN_CPM + 1)).sum(axis=1) >= COEXPR_MIN_CLUSTERS
    lc = lc.loc[expressed]
    annot = prof.loc[keep, ["n_cells", "dominant_cell_class", "purity_cell_class",
                            "dominant_region", "age_pcw_median", "purity_donor"]]
    return lc, annot


COEXPR_MIN_CLUSTERS_PER_CLASS = 5


def context_values(lc: pd.DataFrame, annot: pd.DataFrame, context: str) -> pd.DataFrame | None:
    """The matrix a co-expression context correlates (genes x clusters).

    across_clusters: log2 CPM as is. within_class: each gene's mean per dominant
    cell class subtracted, keeping classes with >= COEXPR_MIN_CLUSTERS_PER_CLASS
    clusters (None if fewer than two such classes' worth of clusters remain).
    """
    if context == "across_clusters":
        return lc
    cls = annot.loc[lc.columns, "dominant_cell_class"].astype(str).to_numpy()
    counts = pd.Series(cls).value_counts()
    keep_cls = set(counts[counts >= COEXPR_MIN_CLUSTERS_PER_CLASS].index)
    cols = np.array([c in keep_cls for c in cls])
    if cols.sum() < 2 * COEXPR_MIN_CLUSTERS_PER_CLASS:
        return None
    Xw = lc.loc[:, cols].to_numpy(float).copy()
    for c in keep_cls:
        m = cls[cols] == c
        Xw[:, m] -= Xw[:, m].mean(axis=1, keepdims=True)
    return pd.DataFrame(Xw, index=lc.index, columns=lc.columns[cols])


def coexpr_contexts(lc: pd.DataFrame, annot: pd.DataFrame) -> dict:
    """{context: (unit-rank rows, matching bins, cluster columns, z-scored values)}.

    across_clusters -- the log2 CPM profiles as they are: co-expression between
                       cell types and states (identity programmes).
    within_class    -- each gene's mean per dominant cell class removed first
                       (classes with >= COEXPR_MIN_CLUSTERS_PER_CLASS clusters):
                       co-variation inside cell types, beyond shared identity.
    Bins (mean level x spread) always come from the raw log2 CPM, so a gene's
    matching bin does not depend on the context.
    """
    out = {}
    for context in ("across_clusters", "within_class"):
        V = context_values(lc, annot, context)
        if V is None:
            continue
        X = V.to_numpy(float)
        out[context] = (unit_rank_rows(X), expression_bins(lc[V.columns]),
                        list(V.columns), zscore_rows(X))
    return out


def unit_rank_rows(X: np.ndarray) -> np.ndarray:
    """Rank each row, centre it and scale to unit length.

    Dot products of the result are Spearman correlations, and the mean pairwise
    correlation of a set follows from the length of its row sum.
    """
    R = rank_rows(np.asarray(X, dtype=float))
    R = R - R.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(R, axis=1, keepdims=True)
    return np.divide(R, nrm, out=np.zeros_like(R), where=nrm > 0)


def expression_bins(lc: pd.DataFrame, n: int = 5) -> np.ndarray:
    """Bin id per gene from quantiles of mean level x spread (n x n bins).

    Highly expressed and highly variable genes correlate more strongly with
    anything, so a fair null for a gene set draws random genes from the same
    level-by-spread bins as its members.
    """
    X = lc.to_numpy(float)
    m, s = X.mean(axis=1), X.std(axis=1)
    qm = np.digitize(m, np.quantile(m, np.linspace(0, 1, n + 1)[1:-1]))
    qs = np.digitize(s, np.quantile(s, np.linspace(0, 1, n + 1)[1:-1]))
    return qm * n + qs


def matched_sets(bins: np.ndarray, idx: np.ndarray, n_sets: int,
                 rng: np.random.Generator) -> np.ndarray:
    """(n_sets x len(idx)) random gene indices, each set matching idx bin for bin.

    Drawn without replacement within a set, so a random set never repeats a
    gene (a repeat would add a spurious r = 1 pair to a coherence null).
    """
    need = pd.Series(bins[idx]).value_counts()
    pools = {b: np.nonzero(bins == b)[0] for b in need.index}
    out = np.empty((n_sets, idx.size), dtype=np.int64)
    for i in range(n_sets):
        col = 0
        for b, k in need.items():
            pool = pools[b]
            take = rng.choice(pool, size=min(k, pool.size), replace=False)
            if take.size < k:          # tiny bin: top up with replacement
                take = np.concatenate([take, rng.choice(pool, size=k - take.size)])
            out[i, col:col + k] = take
            col += k
    return out


def set_sums(Zn: np.ndarray, sets: np.ndarray) -> np.ndarray:
    """Row sums of Zn over each set: (n_sets x k) -> (n_sets x columns)."""
    return set_mean_rows(Zn, sets) * sets.shape[1]


def coherence_from_sums(sums: np.ndarray, k: int) -> np.ndarray:
    """Mean pairwise correlation of unit rows from their sum: (|S|^2 - k) / (k(k-1))."""
    return ((sums ** 2).sum(axis=-1) - k) / (k * (k - 1))


def fmt_p(p: float) -> str:
    return "NA" if not np.isfinite(p) else (f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}")


# ---------------------------------------------------------------------------
# per-stratum test -> v2 x v3 combination (shared by 11-20)
# ---------------------------------------------------------------------------
def level_bins(level, n_bins: int = 10) -> np.ndarray:
    """Quantile bin (0 .. n_bins-1) of each gene's expression level."""
    level = np.asarray(level, dtype=float)
    return np.digitize(level, np.quantile(level, np.linspace(0, 1, n_bins + 1)[1:-1]))


def bin_pools(bins: np.ndarray) -> dict[int, np.ndarray]:
    return {int(b): np.nonzero(bins == b)[0] for b in np.unique(bins)}


def matched_draws(bins: np.ndarray, idx: np.ndarray, n_sets: int,
                  rng: np.random.Generator) -> np.ndarray:
    """(n_sets x len(idx)) random genes, each from its member's bin (with replacement)."""
    pools = bin_pools(bins)
    return np.stack([rng.choice(pools[int(bins[i])], size=n_sets) for i in idx], axis=1)


def null_effect(obs: float, null) -> tuple[float, float, float, float]:
    """(effect in null SDs, two-sided empirical p, null mean, null SD)."""
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size < 2:
        return np.nan, np.nan, np.nan, np.nan
    mu, sd = float(null.mean()), float(null.std())
    p = (np.sum(np.abs(null - mu) >= abs(obs - mu) - 1e-12) + 1) / (null.size + 1)
    return ((obs - mu) / sd if sd > 0 else np.nan), float(p), mu, sd


def combine_chemistries(per: pd.DataFrame, keys: list[str], effect: str = "effect_vs_null_sd",
                        p: str = "perm_p", weight: str | None = None,
                        labels: tuple[str, str] = ("higher", "lower"),
                        carry: tuple[str, ...] = ()) -> pd.DataFrame:
    """v2 x v3 per key: signed Stouffer, BH per dataset (if a key), tier, direction.

    per holds one row per key x chemistry with an effect (signed) and a
    two-sided p; weight names a column whose square root weights each
    chemistry (e.g. the number of age points). Columns in carry are kept
    as <col>_v2 / <col>_v3.
    """
    rows = []
    for k, g in per.groupby(keys, sort=False):
        g = g.drop_duplicates("chemistry").set_index("chemistry")
        if not set(CHEMISTRIES) <= set(g.index):
            continue
        a, b = g.loc["v2"], g.loc["v3"]
        w = [np.sqrt(float(a[weight])), np.sqrt(float(b[weight]))] if weight else [1.0, 1.0]
        Z, pc = signed_stouffer([np.array([a[effect]], float), np.array([b[effect]], float)],
                                [np.array([a[p]], float), np.array([b[p]], float)], w)
        row = dict(zip(keys, k if isinstance(k, tuple) else (k,)))
        for col in carry:
            row[f"{col}_v2"], row[f"{col}_v3"] = a[col], b[col]
        row.update({"effect_v2": a[effect], "p_v2": a[p], "effect_v3": b[effect], "p_v3": b[p],
                    "stouffer_z": float(Z[0]), "combined_p": float(pc[0])})
        rows.append(row)
    comb = pd.DataFrame(rows)
    if comb.empty:
        return comb
    comb["combined_q"] = np.nan
    groups = comb.groupby("dataset").groups.items() if "dataset" in keys else [(None, comb.index)]
    for _, ix in groups:
        comb.loc[ix, "combined_q"] = bh(comb.loc[ix, "combined_p"])
    comb["tier"] = replication_tier(comb.effect_v2, comb.p_v2, comb.effect_v3, comb.p_v3,
                                    comb.combined_q)
    comb["direction"] = np.where(comb.stouffer_z > 0, labels[0], labels[1])
    return comb.sort_values("combined_p", kind="stable").reset_index(drop=True)


def partial_spearman(X: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Spearman of each row of X with y, controlling for z (ranks residualised on z)."""
    R = rank_rows(np.asarray(X, dtype=float))
    ry, rz = stats.rankdata(y), stats.rankdata(z)
    R = R - R.mean(axis=1, keepdims=True)
    ry, rz = ry - ry.mean(), rz - rz.mean()
    zz = float(dot(rz, rz))
    Rres = R - np.outer(dot(R, rz) / zz, rz)
    yres = ry - (float(dot(ry, rz)) / zz) * rz
    den = np.linalg.norm(Rres, axis=1) * np.linalg.norm(yres)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, dot(Rres, yres) / den, np.nan)


def gene_length(dataset: str) -> pd.Series:
    """Genomic span (End - Start, bp) per symbol, from the cross-dataset map."""
    g = csv("_cross_dataset", "gene_id_map.csv", low_memory=False)
    g = g[g["dataset"] == dataset].drop_duplicates("symbol")
    span = pd.to_numeric(g["End"], errors="coerce") - pd.to_numeric(g["Start"], errors="coerce")
    return pd.Series(span.to_numpy(), index=g["symbol"].astype(str).to_numpy()).where(lambda s: s > 0)


def top_markers(namespace: str, grouping: str, n: int = 5) -> dict[str, str]:
    """{group: 'GENE1, GENE2, ...'} from 10_markers (empty if the table is absent)."""
    p = EXPORTS / namespace / "10_markers" / f"top_markers_{grouping}.csv"
    if not p.exists():
        return {}
    t = pd.read_csv(p)
    t["top_group"] = t["top_group"].astype(str)
    t = t.sort_values(["top_group", "rank_in_group"])
    return {g: ", ".join(h["gene"].astype(str).head(n)) for g, h in t.groupby("top_group")}


# Some clusterings were computed per chemistry, reusing labels for different
# cells (cortex Clusters / ClustersModularity / ClustersSurprise); others were
# computed once (cortex leiden_scVI / louvain, human_dev cluster_id). Only the
# latter allow "the same cluster" to be compared between v2 and v3.
SHARED_LABEL_IDENTITY = 0.7


def cluster_label_identity(dataset: str, clustering: str, min_cells: int = COEXPR_MIN_CELLS,
                           n_genes: int = 2000) -> tuple[float, int]:
    """(share of labels whose v2 profile best matches the same label in v3, labels compared).

    Profiles: log2 TMM-CPM of clusters with >= min_cells in both chemistries,
    each gene centred over those clusters, the n_genes most variable genes,
    Pearson between clusters.
    """
    mats = [tmm_log_cpm(group_matrix(ns(dataset, c), f"cluster_{clustering}", "pseudobulk_counts",
                                     min_cells=min_cells)) for c in CHEMISTRIES]
    a, b = mats
    shared = a.columns.intersection(b.columns)
    if len(shared) < 3:
        return 0.0, len(shared)
    genes = a.index.intersection(b.index)
    za, zb = a.loc[genes, shared], b.loc[genes, shared]
    za, zb = za.sub(za.mean(axis=1), axis=0), zb.sub(zb.mean(axis=1), axis=0)
    top = (za.std(axis=1) + zb.std(axis=1)).nlargest(n_genes).index
    R = corr_rows(za.loc[top].T.to_numpy(float), zb.loc[top].T.to_numpy(float))
    return float((R.argmax(axis=1) == np.arange(len(shared))).mean()), len(shared)


# ---------------------------------------------------------------------------
# shared by the part-B analyses (21-28)
# ---------------------------------------------------------------------------
def gene_age_trends(lc: pd.DataFrame, ages) -> pd.DataFrame:
    """Per gene: Spearman with age and its permutation p (exact for few points)."""
    ages = np.asarray(ages, dtype=float)
    rho = spearman_rows(lc.to_numpy(float), ages)
    p, _ = spearman_perm_p(rho, ages)
    return pd.DataFrame({"gene": lc.index.to_numpy(), "rho": rho, "perm_p": p,
                         "mean_log2cpm": lc.mean(axis=1).to_numpy()})


def set_shift_test(values: pd.Series, level: pd.Series, sets: dict[str, list[str]],
                   rng: np.random.Generator, n_random: int = 2000, min_genes: int = 5,
                   n_bins: int = 10) -> pd.DataFrame:
    """Mean of a per-gene statistic over each set, vs random sets matched on expression decile.

    values and level are indexed by gene; genes with a missing value are left out.
    Returns gene_set, n_genes, mean, null_mean, null_sd, effect_vs_null_sd, perm_p.
    """
    ok = values.notna() & level.reindex(values.index).notna()
    v = values[ok].to_numpy(float)
    genes = values.index[ok]
    bins = level_bins(level.reindex(genes).to_numpy(float), n_bins)
    pos = {g: i for i, g in enumerate(genes)}
    rows = []
    for name, members in sets.items():
        idx = np.array(sorted({pos[g] for g in members if g in pos}))
        if idx.size < min_genes:
            continue
        obs = float(v[idx].mean())
        null = set_mean_rows(v[:, None], matched_draws(bins, idx, n_random, rng))[:, 0]
        eff, p, mu, sd = null_effect(obs, null)
        rows.append({"gene_set": name, "n_genes": int(idx.size), "mean": obs, "null_mean": mu,
                     "null_sd": sd, "effect_vs_null_sd": eff, "perm_p": p})
    return pd.DataFrame(rows)


def parse_group(label: str, n: int) -> list[str]:
    """'Radial glia | Telencephalon | 8.0' -> ['Radial glia', 'Telencephalon', '8.0'] (n parts)."""
    parts = [x.strip() for x in str(label).split("|")]
    return parts if len(parts) == n else [str(label)] + [""] * (n - 1)


# ---------------------------------------------------------------------------
# external annotations (part C, analyses 29-34); fetched by
# running_scripts/fetch_annotations.sh into config.ANNOTATIONS_DIR
# ---------------------------------------------------------------------------
def annotations_dir() -> Path:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    import config
    return Path(config.ANNOTATIONS_DIR)


def annotation(rel: str) -> Path | None:
    """Path of an annotation file, or None when it has not been fetched."""
    p = annotations_dir() / rel
    return p if p.exists() and p.stat().st_size > 0 else None


def read_gmt(path: Path) -> dict[str, list[str]]:
    """GMT: name <tab> description <tab> genes... (Enrichr leaves the description empty)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            genes = [g.split(",")[0].strip() for g in parts[2:] if g.strip()]
            if genes:
                out[parts[0].strip()] = list(dict.fromkeys(genes))
    return out


_HGNC: dict[str, dict] = {}


def hgnc_to_ensembl() -> dict[str, str]:
    """UPPER-CASE symbol -> Ensembl id, from HGNC: current symbols, then unambiguous
    previous symbols, then unambiguous aliases (never overriding a current symbol)."""
    if "map" in _HGNC:
        return _HGNC["map"]
    p = annotation("hgnc/hgnc_complete_set.txt")
    out: dict[str, str] = {}
    if p is not None:
        h = pd.read_csv(p, sep="\t", low_memory=False, usecols=["symbol", "prev_symbol", "alias_symbol",
                                                                 "ensembl_gene_id"])
        h = h.dropna(subset=["ensembl_gene_id"])
        for sym, ens in zip(h["symbol"].astype(str), h["ensembl_gene_id"].astype(str)):
            out[sym.upper()] = ens
        for col in ("prev_symbol", "alias_symbol"):
            pairs = [(s.strip().upper(), e) for v, e in zip(h[col], h["ensembl_gene_id"]) if isinstance(v, str)
                     for s in v.split("|") if s.strip()]
            counts = pd.Series([s for s, _ in pairs]).value_counts()
            for s, e in pairs:
                if counts[s] == 1 and s not in out:
                    out[s] = e
    _HGNC["map"] = out
    return out


_SYMBOL_MAP: dict[str, dict[str, str]] = {}


def to_dataset_symbols(genes, dataset: str) -> dict[str, str]:
    """{input symbol: this dataset's symbol} for annotation genes (map_genes, cached)."""
    cache = _SYMBOL_MAP.setdefault(dataset, {})
    todo = [g for g in dict.fromkeys(map(str, genes)) if g not in cache]
    if todo:
        m = map_genes(todo, dataset)
        for inp, sym, how in zip(m["input"], m["symbol"], m["match"]):
            cache[inp] = sym if how != "missing" else ""
    return {g: cache.get(str(g), "") for g in genes if cache.get(str(g), "")}


def gmt_libraries() -> dict[str, dict[str, list[str]]]:
    """{library: {term: genes}} for every .gmt in <annotations>/gmt/."""
    d = annotations_dir() / "gmt"
    return {p.stem: read_gmt(p) for p in sorted(d.glob("*.gmt"))} if d.is_dir() else {}


def tf_list() -> set[str]:
    p = annotation("tf/TF_names_v_1.01.txt")
    return {x.strip() for x in p.read_text().splitlines() if x.strip()} if p else set()


def collectri() -> pd.DataFrame:
    """TF -> target edges (CollecTRI via OmniPath): tf, target, sign (+1 / -1 / 0)."""
    p = annotation("tf/collectri.tsv")
    if p is None:
        return pd.DataFrame(columns=["tf", "target", "sign"])
    t = pd.read_csv(p, sep="\t")
    t = t[~t.source_genesymbol.astype(str).str.contains("COMPLEX|_", regex=True)]
    sign = np.where(t.consensus_stimulation & ~t.consensus_inhibition, 1,
                    np.where(t.consensus_inhibition & ~t.consensus_stimulation, -1, 0))
    return pd.DataFrame({"tf": t.source_genesymbol.astype(str), "target": t.target_genesymbol.astype(str),
                         "sign": sign}).drop_duplicates(["tf", "target"])


def gnomad_constraint() -> pd.DataFrame:
    """Per gene symbol: loeuf (LoF o/e upper CI), pli, mis_z -- MANE / canonical transcript."""
    p = annotation("constraint/gnomad.v4.1.constraint_metrics.tsv")
    if p is None:
        return pd.DataFrame(columns=["loeuf", "pli", "mis_z"])
    c = pd.read_csv(p, sep="\t", usecols=["gene", "mane_select", "canonical", "lof.oe_ci.upper", "lof.pLI",
                                           "mis.z_score"], low_memory=False)
    c["_rank"] = c["mane_select"].astype(str).eq("true") * 2 + c["canonical"].astype(str).eq("true")
    c = c.sort_values("_rank", ascending=False).drop_duplicates("gene")
    out = pd.DataFrame({"loeuf": pd.to_numeric(c["lof.oe_ci.upper"], errors="coerce").to_numpy(),
                        "pli": pd.to_numeric(c["lof.pLI"], errors="coerce").to_numpy(),
                        "mis_z": pd.to_numeric(c["mis.z_score"], errors="coerce").to_numpy()},
                       index=c["gene"].astype(str).to_numpy())
    return out[out.loeuf.notna()]


def disease_genes() -> pd.DataFrame:
    """HPO gene -> disease (OMIM / ORPHA ids): gene, disease_id, association_type."""
    p = annotation("disease/genes_to_disease.txt")
    if p is None:
        return pd.DataFrame(columns=["gene", "disease_id", "association_type"])
    h = pd.read_csv(p, sep="\t")
    return h.rename(columns={"gene_symbol": "gene"})[["gene", "disease_id", "association_type"]]


def ligand_receptor_pairs(min_resources: int = 4) -> pd.DataFrame:
    """Ligand -> receptor pairs: OmniPath interactions whose source is annotated as a ligand
    (secreted or membrane) and target as a receptor (plasma membrane) by >= min_resources
    intercell databases each."""
    pi, pr = annotation("lr/omnipath_interactions.tsv"), annotation("lr/intercell_ligand_receptor.tsv")
    if pi is None or pr is None:
        return pd.DataFrame(columns=["ligand", "receptor"])
    ic = pd.read_csv(pr, sep="\t", low_memory=False)
    ic = ic[ic.entity_type == "protein"]
    # receptors must sit in the plasma membrane; ligands must be secreted or on the membrane
    flags = ic.groupby("genesymbol")[["secreted", "plasma_membrane_transmembrane",
                                      "plasma_membrane_peripheral"]].any()
    membrane = set(flags.index[flags.plasma_membrane_transmembrane | flags.plasma_membrane_peripheral])
    outside = membrane | set(flags.index[flags.secreted])
    ic = ic[((ic.category == "receptor") & ic.genesymbol.isin(membrane))
            | ((ic.category == "ligand") & ic.genesymbol.isin(outside))]
    n = ic.groupby(["category", "genesymbol"]).database.nunique()
    lig = set(n.loc["ligand"][n.loc["ligand"] >= min_resources].index) if "ligand" in n.index.get_level_values(0) else set()
    rec = set(n.loc["receptor"][n.loc["receptor"] >= min_resources].index) if "receptor" in n.index.get_level_values(0) else set()
    it = pd.read_csv(pi, sep="\t")
    it = it[it.source_genesymbol.isin(lig) & it.target_genesymbol.isin(rec)]
    return (pd.DataFrame({"ligand": it.source_genesymbol.astype(str), "receptor": it.target_genesymbol.astype(str)})
            .drop_duplicates().reset_index(drop=True))
