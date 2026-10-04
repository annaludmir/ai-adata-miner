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


def rank_rows(X: np.ndarray) -> np.ndarray:
    return np.apply_along_axis(stats.rankdata, 1, X)


def spearman_rows(X: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Spearman rho of every row of X against the vector x (ties averaged)."""
    R = rank_rows(np.asarray(X, dtype=float))
    rx = stats.rankdata(x)
    R = R - R.mean(axis=1, keepdims=True)
    rx = rx - rx.mean()
    den = np.sqrt((R ** 2).sum(axis=1) * (rx ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, R @ rx / den, np.nan)


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
        null = rxp @ base / np.sqrt((base ** 2).sum() * (rxp ** 2).sum(axis=1))
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
        null = np.abs(np.where(den > 0, R @ rxp.T / den, 0.0))
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


def fmt_p(p: float) -> str:
    return "NA" if not np.isfinite(p) else (f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}")
