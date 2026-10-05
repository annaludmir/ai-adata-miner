# ai-adata-miner

Turn two large scRNA-seq atlases of the developing human brain into a CSV layer
that downstream (and AI-generated) analyses can work on without ever touching a
multi-gigabyte `.h5ad` again.

```
schemas/            stored + regenerated AnnData schemas
config.py           dataset registry and semantic column roles
lib/                io (h5ad -> frames), chunked aggregation, statistics, gene panels
scripts_generated/  00-18, the extraction pipeline (steps 1 & 2)
csv_exports/        CSV outputs, one folder per dataset, each with _manifest.csv
panels/             drop real SFARI / DDG2P exports here to override the seed lists
prompts/            system prompts driving step 3
downstream_analyses/  step 3: analyses over the CSVs, results/ and REPORT.md
running_scripts/    Slurm job scripts for powerslurm
docs/               ANALYSIS_CATALOG.md -- what is extracted and what is possible
```

## Datasets

| key | file | cells × genes |
|---|---|---|
| `cortex` | `Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad` | 297,927 × 33,538 |
| `human_dev` | `human_dev_without_week_5.h5ad` | 1,665,937 × 59,459 |

Paths resolve through `config.py`; override with `AI_ADATA_DATA_ROOT` and
`AI_ADATA_OUT_ROOT`.

Both datasets mix 10x **v2 and v3 chemistry**, which in `human_dev` is nearly
confounded with age (v2 ~6–10 pcw, v3 ~5–5.5 and 11.5–14). Analyses therefore
run **per chemistry by default**, into `csv_exports/<dataset>__v2/` and
`__v3/`; `--chemistry all` pools them instead. Script 18 measures what the
stratification costs — read it before interpreting anything.

## Excluding cells

`exclusions.csv` lists cells to leave out of every analysis, one rule per row:

```
dataset,role,value,reason
human_dev,age,5.0,"File is human_dev_without_week_5 but still holds ... cells at 5.0 pcw"
```

`role` is a semantic role (`age`, `donor`, `sample`, `region`, ...) or a raw
`.obs` column; `dataset` may be `*`. Ages match numerically, and donors match
across the two files' ID styles. Every extraction step applies the rules, and
script 01 records what each one removed in `01_overview/exclusions_applied.csv`.
Step 3 re-applies the age, donor and sample rules to older exports. Use
`--exclusions other.csv` (or `AIM_EXCLUSIONS=...` on the cluster) to swap the
file, and `--exclusions none` to switch exclusions off.

## Gene lists

Put gene lists in one folder, one list per file: a CSV with a `gene` column (or
the list in its first column), or plain text with one gene per line. The file
name, without extension, becomes the list's name. Symbols match
case-insensitively, and Ensembl ids are accepted.

- On the cluster, the default folder is
  `/miridan-data/annaludmir/ndd_gene_modules/data/genes/final_genes_to_run_on`.
- Locally, it is `gene_lists/` in the repo. Git ignores that folder, so the lists
  are never uploaded.
- Set `AIM_GENE_LISTS=/path` to use any other folder.

The lists join the gene panels as group `user_lists`. Stage 2 then always
exports their genes, however lowly expressed, and scripts 11 and 17 score them.
Step-3 analyses 06-08 study them as groups.

GWAS-derived lists name every gene near an associated variant, and neighbouring
genes are often co-regulated. In one bipolar list, 16 histones from the 6p22
cluster made the whole list look like a cell-division program. So in step 3,
lists whose name contains `GWAS` keep one gene per locus: genes on one
chromosome within 1 Mb of each other chain into a single locus, and its most
highly expressed gene is kept. Position in the file is no guide here, because
GWAS lists are usually sorted by coordinate. `AIM_COLLAPSE_LISTS=<regex>` changes which lists this
applies to, and `none` turns it off. 06 reports every list's multi-gene loci.

Step 3 matches list genes in four ways: exact symbol, case-insensitive symbol,
Ensembl id, and finally the other file's annotation. The fourth catches genes
renamed between the two files' annotations; for example, cortex calls a
histone HIST1H1C where human_dev calls it H1-2.

## Running

```bash
pip install -r requirements.txt

# smoke test first -- same code path, 20k cells, a couple of minutes
./run_all.sh --limit-cells 20000

# then the real thing
./run_all.sh
```

### On powerslurm

```bash
cd /miridan-data/annaludmir/ai-adata-miner/running_scripts
LIMIT_CELLS=20000 sbatch --export=ALL slurm_full_pipeline.sh   # smoke test first
./submit_all.sh                                                # then the full chain
```

`submit_all.sh` runs the three stages as a dependency chain, with the two
datasets' pseudobulk passes concurrent. See
[`running_scripts/README.md`](running_scripts/README.md) for the resource
layout, the environment knobs, and the cluster gotchas.

### Individual steps

Individual steps take the same flags:

```bash
python3 scripts_generated/09_pseudobulk.py --dataset human_dev --chunk-size 100000
```

Scripts 00–08 and 13 read only `.obs` / `.var` / `.obsm` / `.varm` and finish in
minutes. Script 09 is the long one: it streams `X` once and fills every grouping's
aggregator in that single pass. Scripts 10 and 11 then work purely on 09's CSVs.

## Step 3: downstream analyses

Runs on CSVs only, so on a laptop, in under a minute:

```bash
./downstream_analyses/run_all.sh          # PYTHON=... to pick an interpreter
```

Findings land in [`downstream_analyses/REPORT.md`](downstream_analyses/REPORT.md),
with per-analysis method, limitations and tables under `downstream_analyses/results/`.
Start with `01_data_audit`: it sets out what the data can support (donors as the
replicate unit, v2/v3 as independent donor sets, cortex as a subset of human_dev).
[`downstream_analyses/README.md`](downstream_analyses/README.md) explains each
analysis and how to add one.

## Design notes

- **No scanpy, and no full `AnnData` is ever loaded.** `lib/io_utils.XReader`
  slices CSR row-blocks straight out of HDF5, so streaming 1.67M cells never
  materialises `obsm` (`Factors` alone would be gigabytes).
- **No column name is hard-coded.** Scripts ask for a semantic *role*
  (`cell_class`, `donor`, `cluster`…) and `config.COLUMN_ROLES` resolves it to
  whichever column that file actually has. This is what lets one script run on
  two differently-annotated datasets, and why a missing column skips an analysis
  with a logged note instead of crashing the run.
- **Every CSV is registered** in `csv_exports/<dataset>/_manifest.csv` with its
  shape, columns and a description, so step 3 discovers inputs instead of
  guessing filenames.
- **The file is the source of truth, not the schema.** Script 00 regenerates the
  schema from the data and warns when it disagrees with `schemas/*.json`.

See [`docs/ANALYSIS_CATALOG.md`](docs/ANALYSIS_CATALOG.md) for the full list of
what is extracted, what is queued for step 3, and the statistical cautions that
apply to all of it.
