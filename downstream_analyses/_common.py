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
                         label: str = "panel") -> pd.DataFrame:
    """Does each gene set score higher in one cell class than in the others?

    Pseudobulk counts per (cell class, age point) -> log2 TMM-CPM; genes >= 5
    CPM on average; each gene Z-scored across all columns; a set's score in a
    column is the mean Z of its genes. For class c, at every age point, the
    difference between c and the mean of the other classes is taken, and T is
    its mean over age points (age points are donors). The null is T for random
    sets drawing one gene from each member's expression decile.
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
    members_of_bin = [np.nonzero(bins == b)[0] for b in range(n_bins)]

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
    and finally the OTHER file's annotation -- a symbol known there is carried
    over through the shared Ensembl id. The two files use different symbol
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
