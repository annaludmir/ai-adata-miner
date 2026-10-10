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
stage 1  metadata, both datasets       16G   1h   scripts 00-08, 13, 18, 20, 21
   |
   +--> stage 2  pseudobulk cortex v2    128G  6h  scripts 09,12,19,22,23─┐
   +--> stage 2  pseudobulk cortex v3    128G  6h  scripts 09,12,19,22,23─┤ all
   +--> stage 2  pseudobulk human_dev v2 128G  6h  scripts 09,12,19,22,23─┤ concurrent
   +--> stage 2  pseudobulk human_dev v3 128G  6h  scripts 09,12,19,22,23─┘
   |
   v  (afterok on all four)
stage 3  derived analyses              16G   3h   scripts 10,11,14-17
   |
   v  (afterok)
stage 4  step-3 downstream analyses     8G   2h   downstream_analyses/ + REPORT.md
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

## External annotations (step-3 part C)

Analyses 29-34 read GO / Reactome / KEGG gene sets, TF lists and targets, gnomAD
constraint, HPO disease genes, ligand-receptor pairs, HGNC symbols and BrainSpan.
Download them once, on the login node (light, ~0.4 GB, about two minutes):

```bash
./fetch_annotations.sh            # -> /miridan-data/annaludmir/aim_annotations, with MANIFEST.tsv
```

Without them those analyses say so and skip. Once the HGNC table is there,
gene-list matching everywhere also uses HGNC previous symbols.

## Spectra gene programmes (scripts 24-25, analysis 36)

Spectra fits gene programmes guided by prior gene sets (lists, seed panels, core
GO processes, class markers). It runs separately from the four stages, one
dataset x chemistry per job, after stages 1-2 and `fetch_annotations.sh`.

Install once into the env (the `setuptools` pin is needed: Spectra still imports
`pkg_resources`):

```bash
/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new/bin/python -m pip install scSpectra "setuptools<81"
```

Call the env's python by path: `$AIM_ENV` is set only inside the Slurm scripts,
and `conda activate` on the login node leaves the shared mamba python and pip
first on PATH. That env's python also reads `~/.local/lib/python3.10`, ahead of
its own packages, so anything installed with `pip --user` changes every stage.
Check what a job will import with
`$ENV/bin/python -c "import anndata; print(anndata.__version__, anndata.__file__)"`.

CPU (the published model; 2-6 h per stratum on 16 cores):

```bash
for ds in cortex human_dev; do for ch in v2 v3; do
  sbatch --export=ALL,DATASET=$ds,CHEMISTRY=$ch slurm_spectra.sh
done; done
```

GPU (the package's minibatched `Spectra_gpu`, which its authors mark as in
development; needs a CUDA build of torch in the env):

```bash
sbatch --partition=<gpu partition> --gres=gpu:1 \
       --export=ALL,DATASET=human_dev,CHEMISTRY=v2,SPECTRA_BACKEND=gpu slurm_spectra.sh
```

The subsample goes to `AIM_WORK` (default `/miridan-data/annaludmir/aim_work`),
the results to `csv_exports/<ns>/25_spectra/`. `SKIP_INPUT=true` refits without
redrawing the subsample; `SPECTRA_EPOCHS` overrides the epoch count.

Spectra's own training rule lowers the learning rate after every 3 epochs whose
loss did not fall and stops after 18, counted over the whole run; the noisy
early epochs fill that count, so one fit ended at epoch 431 and the others spent
most of their 5000 epochs at a rate of 0.001. Script 25 replaces it with a
plateau rule: the same rates, each lowered only after `SPECTRA_PLATEAU` (default
50) epochs without a new best loss, stopping when the last rate stops improving
or at `SPECTRA_EPOCHS` (default 10000). `SPECTRA_PLATEAU=0` restores Spectra's
rule. The job log and `25_spectra/fit_summary.csv` give the epochs run, why
training stopped and the rate steps; `training_trace.csv` the loss per epoch.

**Second fit, core prior.** With the gene lists in the prior they dominate the
global gene-gene graph, so even the free factors are pulled into them and the
lists cannot be tested against the programmes. `SPECTRA_PRIOR=core` fits GO
processes, cell-cycle panels and class markers only, with 15 free global
factors, into `25_spectra_core/` (work files `*_core`); analysis 36 tests the
NDD lists against those programmes. Reuse the subsamples:

```bash
for ds in cortex human_dev; do for ch in v2 v3; do
  sbatch --export=ALL,DATASET=$ds,CHEMISTRY=$ch,SKIP_INPUT=true,SPECTRA_PRIOR=core slurm_spectra.sh
done; done
```

To rewrite the per-group scores of finished fits (e.g. after a grouping is
added) without refitting, run on a compute node:

```bash
srun --account=miridan-users_v2 --partition=power-general-public-pool --qos=public --mem=16G --time=00:30:00 \
  /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new/bin/python scripts_generated/25_spectra_fit.py --rescore
```

(add `--prior core` for the core-prior fits)
Fit all four strata with the same settings. Rerun
step 3 (or just analysis 36) afterwards.

## Knobs

All read from the environment, all forwarded through the chain by `submit_all.sh`.

| variable | default | notes |
|---|---|---|
| `DATASET` / `DATASETS` | `all` / `cortex human_dev` | which dataset(s) |
| `CHEMISTRY` | `each` | `each` runs v2 and v3 apart; `all` pools them |
| `LIMIT_CELLS` | unset | first N cells — smoke test; outputs are partial |
| `AIM_EXCLUSIONS` | `exclusions.csv` | cells to leave out (see the top-level README); `none` disables |
| `AIM_GENE_LISTS` | `/miridan-data/annaludmir/ndd_gene_modules/data/genes/final_genes_to_run_on` | folder of gene lists (one per file) analysed in step 3 and always exported by stage 2 |
| `AIM_COLLAPSE_LISTS` | `GWAS` | regex of list names collapsed to one gene per locus in step 3; `none` disables |
| `AIM_DOWNSTREAM_OUT` | `/miridan-data/annaludmir/aim_downstream/results` | where stage 4 writes step-3 results |
| `AIM_ANNOTATIONS` | `/miridan-data/annaludmir/aim_annotations` | external annotation files for step-3 analyses 29-34; fill it once with `./fetch_annotations.sh` on the login node |
| `AIM_STRICT` | `false` | `true` reruns with `exclusions_strict.csv` into `csv_exports_strict/` and `.../aim_downstream_strict/results` |
| `AIM_STRICT_RESULTS` | `.../aim_downstream_strict/results` | strict-QC results that step-3 analysis 28 compares against |
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
