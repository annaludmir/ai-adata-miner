# Running the pipeline on powerslurm

Follows the conventions of `ndd_gene_modules/running_scripts/`: public pool,
`miridan-users_v2`, logs in `/miridan-data/annaludmir/jobs_output/%j.{out,err}`,
the `mamba-1.5.8` + `jupyter-scanpy_new` environment, and the EXIT-trap summary
file from `python_cell_cycle_annotation.sh`.

## Order of operations

```bash
cd /miridan-data/annaludmir/ai-adata-miner/running_scripts
bash check_env.sh                                              # 1. is the env usable?  (~1 min)
LIMIT_CELLS=20000 sbatch --export=ALL slurm_full_pipeline.sh   # 2. smoke test          (~minutes)
./submit_all.sh                                                # 3. the full run        (hours)
```

## First: check the environment

`check_env.sh` confirms the five required packages, resolves `anndata.read_elem`,
imports the repo modules and opens both h5ad files. See
[Creating the environment](#creating-the-environment) — in short, the existing
`jupyter-scanpy_new` almost certainly already works and nothing needs building.

## Then: a smoke test

Same code path, 20k cells, minutes instead of hours. Do this before committing
the full allocation.

```bash
cd /miridan-data/annaludmir/ai-adata-miner/running_scripts
LIMIT_CELLS=20000 sbatch --export=ALL slurm_full_pipeline.sh
```

Then read `/miridan-data/annaludmir/jobs_output/<jobid>_summary.txt`. It lists
every step with its runtime, anything that failed, the CSV count, and the
warnings worth acting on (schema mismatches, confounds, skipped groupings).

## Then: the full run

```bash
./submit_all.sh          # run on the LOGIN node, not via sbatch
```

That submits a four-stage chain:

```
stage 1  metadata, both datasets       16G   1h   scripts 00-08, 13, 18
   |
   +--> stage 2  pseudobulk cortex v2    128G  6h  scripts 09, 12  ─┐
   +--> stage 2  pseudobulk cortex v3    128G  6h  scripts 09, 12  ─┤ all
   +--> stage 2  pseudobulk human_dev v2 128G  6h  scripts 09, 12  ─┤ concurrent
   +--> stage 2  pseudobulk human_dev v3 128G  6h  scripts 09, 12  ─┘
   |
   v  (afterok on all four)
stage 3  derived analyses              16G   1h   scripts 10,11,14-17
   |
   v  (afterok)
stage 4  step-3 downstream analyses     8G  30m   downstream_analyses/ + REPORT.md
```

Stage 4 writes **outside the repo**, to `AIM_DOWNSTREAM_OUT` (default
`/miridan-data/annaludmir/aim_downstream/results`, with `REPORT.md` one level
up). `REPORT.md`, the `SUMMARY.md` files and the figures are tracked in git, so
writing them into the cluster checkout would make the next `git pull` conflict.
Read them there, or copy them back with `rsync`. `DOWNSTREAM=false ./submit_all.sh`
stops after stage 3. To rerun step 3 alone after changing an analysis:

```bash
sbatch slurm_04_downstream.sh
```

Stage 2 fans out over dataset x chemistry because the default `--chemistry each`
needs one pass per chemistry; running them concurrently spends cluster slots
instead of wall time. `CHEMISTRY=all ./submit_all.sh` halves the job count and
pools the chemistries, at the cost of confounding age with the v2->v3 switch.

The split exists for three reasons. Stage 1 never touches the count matrix, so
holding stage-2 memory for it is waste. The two stage-2 jobs read different files and
write different folders, so running them together roughly halves wall time
(human_dev is ~5x the cells of cortex). And stage 3 reads only `csv_exports/`,
so re-running it after changing a panel or a K sweep costs minutes and never
re-reads the matrix.

| script | what it does |
|---|---|
| `check_env.sh` | verifies AIM_ENV and that both h5ad files open — run this first |
| `submit_all.sh` | login-node orchestrator; submits the chain with `afterok` |
| `slurm_01_metadata.sh` | stage 1 — obs/var/obsm only |
| `slurm_02_pseudobulk.sh` | stage 2 — the single streaming pass over X |
| `slurm_03_derived.sh` | stage 3 — CSV-only analyses |
| `slurm_full_pipeline.sh` | all 18 steps in one job (simpler; no concurrency) |
| `_common.sh` | sourced helper: env, per-step runner, summary trap |

## Knobs

All read from the environment, all forwarded through the chain by `submit_all.sh`.

| variable | default | notes |
|---|---|---|
| `DATASET` / `DATASETS` | `all` / `cortex human_dev` | which dataset(s) |
| `CHEMISTRY` | `each` | `each` runs v2 and v3 apart; `all` pools them |
| `LIMIT_CELLS` | unset | first N cells — smoke test; outputs are partial |
| `AIM_EXCLUSIONS` | `exclusions.csv` | cells to leave out (see the top-level README); `none` disables |
| `AIM_DOWNSTREAM_OUT` | `/miridan-data/annaludmir/aim_downstream/results` | where stage 4 writes step-3 results |
| `DOWNSTREAM` | `true` | `false` makes `submit_all.sh` stop after stage 3 |
| `CHUNK_SIZE` | `50000` | cells per streaming chunk (stage 2) |
| `TOP_GENES` | `12000` | genes exported in pseudobulk; panel genes always added |
| `AIM_ROOT` | `/miridan-data/annaludmir/ai-adata-miner` | repo location |
| `AIM_ENV` | `.../conda-envs/jupyter-scanpy_new` | mamba prefix |
| `AI_ADATA_DATA_ROOT` | `.../ndd_gene_modules/data` | where the h5ad files live |
| `AI_ADATA_OUT_ROOT` | `$AIM_ROOT/csv_exports` | where CSVs are written |
| `MAIL_SUMMARY` | `false` | `true` → watch the final job and mail its summary |

```bash
DATASET=human_dev CHUNK_SIZE=25000 sbatch --export=ALL slurm_02_pseudobulk.sh
TOP_GENES=20000 MAIL_SUMMARY=true ./submit_all.sh
```

## Things that will bite you

- **`--export=ALL` does not carry variables through a dependency chain.**
  `submit_all.sh` forwards each one explicitly. If you submit a stage by hand,
  pass `--export=ALL` *and* set the variable in the same command.
- **Compute nodes accept mail but silently fail to relay it.** Same finding as
  `python_cell_cycle_annotation.sh`. `EMAIL_OUTPUT` is therefore `false` and the
  summary *file* is the reliable copy. To get mail, use the login-node helper:
  `MAIL_SUMMARY=true ./submit_all.sh`, which runs `mail_job_summary.sh`.
- **A failed stage leaves its dependents queued as `DependencyNeverSatisfied`.**
  They never run and Slurm does not always clear them — `scancel` them (the
  command is printed at submission).
- **A failed *step* does not abort its stage.** The runner records it and
  carries on, so one broken script does not discard the hours step 09 already
  spent. The stage still exits 1, so `afterok` dependents correctly do not fire.
- **Slurm runs batch scripts from a spool copy**, so `$BASH_SOURCE` is not the
  repo. The scripts resolve `_common.sh` through `$AIM_ROOT` instead; if you
  move the repo, set `AIM_ROOT`.
- **Resource figures are estimates, not measurements.** They were never run
  against the real files. Check `seff <jobid>` after the first real run and
  tune — in particular `--mem` for stage 2, which is the only one that matters.

## Creating the environment

`AIM_ENV` points at a mamba environment prefix. **You most likely do not need to
create one.**

### First, check the environment you already have

`AIM_ENV` defaults to `/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new`,
the env the `ndd_gene_modules` jobs use. The pipeline needs exactly five
packages — numpy, pandas, scipy, h5py, anndata — and that env already has them
(scanpy depends on anndata and h5py). So the first step is to confirm, not to
build:

```bash
cd /miridan-data/annaludmir/ai-adata-miner/running_scripts
bash check_env.sh            # ~1 minute, light enough for the login node
```

It reports the five package versions, resolves `anndata.read_elem` (which has
moved twice across versions — the pipeline handles all three locations), imports
the repo's own modules, and then opens both h5ad files and reads their shape and
`.obs`. It also tells you if the real files disagree with `schemas/`. If it ends
in `READY`, submit jobs — there is nothing to create.

### Only if that fails: build a dedicated env

```bash
module load mamba/mamba-1.5.8
mamba env create -p /miridan-data/annaludmir/conda-envs/ai-adata-miner \
                 -f /miridan-data/annaludmir/ai-adata-miner/environment.yml

AIM_ENV=/miridan-data/annaludmir/conda-envs/ai-adata-miner bash check_env.sh
```

Then pass it to every job:

```bash
export AIM_ENV=/miridan-data/annaludmir/conda-envs/ai-adata-miner
./submit_all.sh                                   # forwards AIM_ENV through the chain
sbatch --export=ALL slurm_02_pseudobulk.sh        # or per job
```

Or edit the default in `_common.sh` once.

This env is deliberately small: no scanpy. The extraction layer reads h5ad
through h5py and `anndata.read_elem` and never constructs an AnnData object, so
scanpy's dependency tree is not needed and the solve is quick.

### Why a separate env might be worth it anyway

The shared `jupyter-scanpy_new` carries pyscenic, pydeseq2, gseapy and a
numpy-compat patch (`fix_pyscenic_numpy_compat.py`). Anything that upgrades
numpy there to satisfy one of those can break the others. A dedicated env with
five pinned packages removes this pipeline from that blast radius. Not urgent —
just the reason you might choose to.
