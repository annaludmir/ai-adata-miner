#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_meta
#SBATCH --mem=16G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Stage 1 of 3 -- everything that never touches the count matrix.
#
# Scripts 00-08 and 13 read only .obs / .var / .obsm / .varm, so they finish in
# minutes even on the 1.67M-cell file and need nowhere near the memory stage 2
# does. Running them first means the inventory, the confound report and the
# factor modules are on disk before the expensive streaming starts -- and if the
# schema disagrees with the file, you find out in minutes rather than hours.
#
# Measured (seff 22526172, both datasets, full data): 2.6 GB peak, 4 min wall,
# one core ~87% busy. 16G / 1h is >6x headroom on memory and ~15x on time.
#
#   sbatch slurm_01_metadata.sh
#   DATASET=cortex sbatch --export=ALL slurm_01_metadata.sh
#   LIMIT_CELLS=20000 sbatch --export=ALL slurm_01_metadata.sh   # smoke test

set -uo pipefail    # not -e: a failed step is recorded and the rest still run

AIM_STAGE="01_metadata"
# Slurm runs this from a spool copy on the compute node, so BASH_SOURCE does not
# point at the repo. Resolve the helper from an absolute root instead.
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup

aim_run 00_inventory_and_schema
aim_run 01_obs_overview_and_qc
aim_run 02_composition
aim_run 03_cellcycle
aim_run 04_cluster_profiles
aim_run 05_covariate_confounds
aim_run 06_embeddings
aim_run 07_factor_modules      # needs varm/Loadings -- human_dev only, skips cortex
aim_run 08_factor_activity
aim_run 18_chemistry_comparability   # pooled by design: what survives stratification

# 13 compares the two datasets, so it only means anything with both of them.
if [[ "$DATASET" == "all" ]]; then
  aim_run 13_cross_dataset_keys
else
  echo ""
  echo "(skipping 13_cross_dataset_keys: needs DATASET=all)"
fi

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
