#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_derived
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=03:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Stage 3 of 3 -- analyses that read only csv_exports/, never the h5ad.
#
# This is the stage worth re-running: changing a panel, a K sweep or a GSEA
# permutation count costs minutes here and does not touch the matrix again.
# It does require stage 2 to have finished, since every script here reads
# 09_pseudobulk/ (17 additionally reads 10_markers/, so order matters).
#
# Measured (seff 22526178, both datasets x both chemistries, 8 CPUs): 1.8 GB
# peak, 8 min wall. Job 22580830 with 2 CPUs and unpinned BLAS threads ran
# ~20x slower and hit a 1h limit inside step 17; _common.sh now pins thread
# counts to SLURM_CPUS_PER_TASK. 3h is margin until the next seff confirms
# the fix -- then it can come back down. Raise --time if you push
# N_PERMUTATIONS in 17 well past 1000.
#
#   sbatch slurm_03_derived.sh
#   sbatch --dependency=afterok:<stage2_jobid> slurm_03_derived.sh

set -uo pipefail

AIM_STAGE="03_derived"
# Slurm runs this from a spool copy on the compute node, so BASH_SOURCE does not
# point at the repo. Resolve the helper from an absolute root instead.
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup

if ! ls "${AI_ADATA_OUT_ROOT}"/*/09_pseudobulk/*__mean_lognorm.csv >/dev/null 2>&1; then
  echo "ERROR: no pseudobulk CSVs under ${AI_ADATA_OUT_ROOT}."
  echo "       Run slurm_02_pseudobulk.sh first."
  exit 1
fi

aim_run 10_marker_specificity
aim_run 11_gene_panels
aim_run 14_normalization_diagnostics
aim_run 15_sample_relationships
aim_run 16_expression_patterns
aim_run 17_gsea_panels          # reads 10_markers/, so it runs after 10

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
