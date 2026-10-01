# Running the pipeline on powerslurm

Follows the conventions of `ndd_gene_modules/running_scripts/`: public pool,
`miridan-users_v2`, logs in `/miridan-data/annaludmir/jobs_output/%j.{out,err}`,
the `mamba-1.5.8` + `jupyter-scanpy_new` environment, and the EXIT-trap summary
file from `python_cell_cycle_annotation.sh`.

## First: a smoke test

Same code path, 20k cells, minutes instead of hours. Do this before committing
a 500G allocation for two days.

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

That submits a three-stage chain:

```
stage 1  metadata, both datasets      250G  12h   scripts 00-08, 13
   |
   +--> stage 2  pseudobulk cortex    500G  48h   scripts 09, 12   ─┐ concurrent
   +--> stage 2  pseudobulk human_dev 500G  48h   scripts 09, 12   ─┘
   |
   v  (afterok on both)
stage 3  derived analyses             200G  12h   scripts 10,11,14-17
```

The split exists for three reasons. Stage 1 never touches the count matrix, so
holding 500G for it is waste. The two stage-2 jobs read different files and
write different folders, so running them together roughly halves wall time
(human_dev is ~5x the cells of cortex). And stage 3 reads only `csv_exports/`,
so re-running it after changing a panel or a K sweep costs minutes and never
re-reads the matrix.

| script | what it does |
|---|---|
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
| `LIMIT_CELLS` | unset | first N cells — smoke test; outputs are partial |
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
