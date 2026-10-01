# ai-adata-miner

Turn two large scRNA-seq atlases of the developing human brain into a CSV layer
that downstream (and AI-generated) analyses can work on without ever touching a
multi-gigabyte `.h5ad` again.

```
schemas/            stored + regenerated AnnData schemas
config.py           dataset registry and semantic column roles
lib/                io (h5ad -> frames), chunked aggregation, statistics, gene panels
scripts_generated/  00-13, the extraction pipeline (steps 1 & 2)
csv_exports/        CSV outputs, one folder per dataset, each with _manifest.csv
panels/             drop real SFARI / DDG2P exports here to override the seed lists
prompts/            system prompts driving step 3
downstream_analyses/  step 3 outputs
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
