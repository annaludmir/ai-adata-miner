"""Central configuration: dataset registry, column semantics, analysis defaults.

The two h5ad files are structurally very different (different gene-id spaces,
different clustering columns, different obsm/layers), so no script hard-codes a
column name.  Each script asks `lib.io_utils.resolve_role` for a semantic role
("cluster", "age", "donor", ...) and takes the first candidate that actually
exists in the file.  That is what lets one script run on both datasets.

Override paths with env vars AI_ADATA_DATA_ROOT / AI_ADATA_OUT_ROOT.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get(
    "AI_ADATA_DATA_ROOT", "/miridan-data/annaludmir/ndd_gene_modules/data"))
CSV_EXPORTS = Path(os.environ.get("AI_ADATA_OUT_ROOT", REPO_ROOT / "csv_exports"))

# Cells left out of every analysis: one rule per row (dataset, role, value,
# reason) in exclusions.csv. Override with AIM_EXCLUSIONS=<path> or the
# --exclusions flag; the value 'none' disables exclusions entirely.
EXCLUSIONS_FILE = os.environ.get("AIM_EXCLUSIONS", str(REPO_ROOT / "exclusions.csv"))

# Folder of user gene lists (one list per file: CSV with a 'gene' column -- or
# its first column -- or plain text, one gene per line; the file stem names the
# list). They join the panels as group 'user_lists', so script 09 always
# exports their genes and 11/17 score them; step 3 analyses them as groups.
# Override with AIM_GENE_LISTS. The repo-local default is git-ignored.
_CLUSTER_GENE_LISTS = Path("/miridan-data/annaludmir/ndd_gene_modules/data/genes/final_genes_to_run_on")
GENE_LISTS_DIR = Path(os.environ.get(
    "AIM_GENE_LISTS",
    _CLUSTER_GENE_LISTS if _CLUSTER_GENE_LISTS.exists() else REPO_ROOT / "gene_lists"))
SCHEMA_DIR = REPO_ROOT / "schemas"

# ---------------------------------------------------------------------------
# Dataset registry
# ---------------------------------------------------------------------------
DATASETS: dict[str, dict] = {
    "cortex": {
        "h5ad": DATA_ROOT / "Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad",
        "label": "Developing human cortex, EMX1 lineage (QC-passed, Post-M)",
        "schema": SCHEMA_DIR / "cortex_schema.json",
        # var index holds gene SYMBOLS here
        "gene_id_space": "symbol",
        "expected_n_obs": 297_927,
        "expected_n_vars": 33_538,
    },
    "human_dev": {
        # The cell-cycle-annotated version: same atlas plus per-cell
        # CellCyclePhase (G1/S/G2M/Post-M/Non-cycling), cycling_score and
        # Age_Chemistry. The plain human_dev_without_week_5.h5ad lacks these.
        "h5ad": DATA_ROOT / "human_dev_without_week_5_cc_annotated.h5ad",
        "label": "Whole developing human brain atlas (cell-cycle annotated)",
        "schema": SCHEMA_DIR / "human_dev_schema.json",
        # var index holds VERSIONED Ensembl accessions here (ENSG...*.N)
        "gene_id_space": "accession_versioned",
        "expected_n_obs": 1_665_937,
        "expected_n_vars": 59_459,
    },
}

# ---------------------------------------------------------------------------
# Semantic roles -> candidate obs columns, in priority order.
# First existing candidate wins.  Unknown roles simply resolve to None and the
# dependent analysis is skipped with a logged note rather than crashing.
# ---------------------------------------------------------------------------
COLUMN_ROLES: dict[str, list[str]] = {
    "cell_class":   ["CellClass", "classes", "cell_type", "CellType"],
    "cell_type_id": ["cell_type_ontology_term_id"],
    # an older or coarser cell-type label kept alongside cell_class (cortex 'classes');
    # 02_composition cross-tabulates it against cell_class for the annotation audit
    "alt_cell_class": ["classes"],
    "age":          ["Age", "age"],
    "age_text":     ["Agetext", "Subset"],
    "region":       ["Region", "Tissue"],
    "subregion":    ["Subregion", "Shortname"],
    "dissection":   ["dissection", "Tissue"],
    "donor":        ["donor_id", "Donor"],
    "sample":       ["sample_id", "SampleID", "Name"],
    "chemistry":    ["Chemistry"],
    "assay":        ["assay_ontology_term_id"],
    "sex":          ["sex_ontology_term_id", "Sex"],
    "dev_stage":    ["development_stage_ontology_term_id"],
    "cyclephase":   ["CellCyclePhase"],
    "age_chem":     ["Age_Chemistry"],
    "transcriptome": ["Transcriptome"],
    "batch":        ["batch"],
}

# Clustering columns are plural: every one present gets profiled separately,
# because comparing alternative partitions of the same cells is itself a result.
CLUSTER_ROLES: list[str] = [
    "cluster_id", "Clusters", "ClustersModularity", "ClustersSurprise",
    "leiden_scVI", "louvain", "louvain3", "leiden",
]

# Numeric QC columns.  Several are stored as *categorical strings* in
# human_dev (total_genes, total_UMIs) -- io_utils coerces them back to numbers.
QC_NUMERIC_ROLES: dict[str, list[str]] = {
    "frac_mito":        ["fraction_mitochondrial", "MT_ratio"],
    "frac_unspliced":   ["fraction_unspliced", "unspliced_ratio"],
    "n_genes":          ["total_genes", "NGenes"],
    "total_umis":       ["total_UMIs", "TotalUMI"],
    "total_rna":        ["TotalRNA"],
    "cell_cycle_score": ["cell_cycle_score", "CellCycle"],
    "cycling_score":    ["cycling_score", "CellCycleFraction"],
    "cc_g1":            ["CellCycle_G1"],
    "cc_s":             ["CellCycle_S"],
    "cc_g2m":           ["CellCycle_G2M"],
    "doublet_score":    ["DoubletFinderScore", "doublet_score", "scrublet_score"],
}

# .var columns carrying gene identity.  NOTE: in human_dev the column literally
# named `gene_symbol` holds Ensembl accessions -- `Gene` is the real symbol
# column.  Order below encodes that correction.
VAR_SYMBOL_CANDIDATES = ["Gene", "gene_name", "symbol", "feature_name", "gene_symbol"]
VAR_ACCESSION_CANDIDATES = ["Accession", "gene_ids", "ensembl_id"]

# Embeddings worth summarising, in preference order.
EMBEDDING_KEYS = ["X_umap", "UMAP", "X_Embedding", "X_scVI", "TSNE", "PCA"]

# Lineage pseudotime (script 20): the latent space each file stores, the
# dorsal neurogenic lineage in order, and -- for the whole-brain atlas -- the
# region whose cells form that lineage (other regions have other lineages).
PSEUDOTIME_EMBEDDING = {"cortex": "X_scVI", "human_dev": "Factors"}
LINEAGE_CLASSES = ["Radial glia", "Neuronal IPC", "Neuroblast", "Neuron"]
PSEUDOTIME_REGIONS = {"human_dev": ["Telencephalon"]}
PSEUDOTIME_PCS = 10            # latent dims are standardised, then reduced to this many PCs
PSEUDOTIME_BINS = 20           # equal-width bins of pseudotime in [0, 1] for pseudobulks

# Label harmonisation so the two datasets can be compared in step 3.
LABEL_SYNONYMS: dict[str, dict[str, str]] = {
    "cyclephase": {"PostM": "Post-M", "Post-M": "Post-M", "post-m": "Post-M"},
}

# ---------------------------------------------------------------------------
# Analysis defaults
# ---------------------------------------------------------------------------
CHUNK_SIZE = 50_000            # cells per pass when streaming X from disk
TARGET_SUM = 1e4               # CP10K normalisation before log1p
MIN_CELLS_PER_GROUP = 20       # thin groups are flagged in the CSV, never dropped
TOP_GENES_PSEUDOBULK = 12_000  # cap on wide gene CSVs; ranked by total UMIs
TOP_GENES_PER_FACTOR = 100     # top loadings exported per factor/module
MAX_GROUPS_WIDE = 1_000        # refuse absurdly wide pseudobulk CSVs
RANDOM_SEED = 0


# Both datasets mix 10x v2 and v3 chemistry, and in human_dev the two cover
# almost disjoint age ranges (v2 ~6-10 pcw, v3 ~5-5.5 and 11.5-14), so a
# "developmental" change measured across pooled cells is partly the chemistry
# switch. Analyses therefore run per chemistry by default; script 05 measures
# the actual overlap on the real file.
CHEMISTRIES = ("v2", "v3")


def namespace(key: str, chemistry: str | None = None) -> str:
    """csv_exports/ folder name for a dataset, optionally chemistry-stratified."""
    return key if chemistry is None else f"{key}__{chemistry}"


def base_dataset(ns: str) -> str:
    """Inverse of namespace(): 'cortex__v2' -> 'cortex'."""
    return ns.split("__")[0]


def dataset(key: str) -> dict:
    key = base_dataset(key)
    if key not in DATASETS:
        raise SystemExit(f"Unknown dataset {key!r}. Known: {', '.join(DATASETS)}")
    return DATASETS[key]


def out_dir(key: str, subdir: str | None = None) -> Path:
    p = CSV_EXPORTS / key if subdir is None else CSV_EXPORTS / key / subdir
    p.mkdir(parents=True, exist_ok=True)
    return p
