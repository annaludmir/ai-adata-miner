#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_full
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# All 18 steps in one job. Simpler than the three-stage chain, and the right
# choice for a smoke test or a rerun you expect to succeed. The trade-off: the
# whole thing holds a 500G allocation for its entire wall time, including the
# cheap metadata steps, and a failure in step 09 leaves you resubmitting
# everything. For a first full run prefer submit_all.sh.
#
#   sbatch slurm_full_pipeline.sh
#   LIMIT_CELLS=20000 sbatch --export=ALL slurm_full_pipeline.sh   # smoke test
#   DATASET=cortex    sbatch --export=ALL slurm_full_pipeline.sh

set -uo pipefail

AIM_STAGE="full_pipeline"
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup

# Order matters: 10/11/14/15/16 read 09's CSVs, and 17 reads 10's.
aim_run 00_inventory_and_schema
aim_run 01_obs_overview_and_qc
aim_run 02_composition
aim_run 03_cellcycle
aim_run 04_cluster_profiles
aim_run 05_covariate_confounds
aim_run 06_embeddings
aim_run 07_factor_modules
aim_run 08_factor_activity
aim_run 09_pseudobulk --top-genes "$TOP_GENES"
aim_run 10_marker_specificity
aim_run 11_gene_panels
aim_run 12_splicing_layers
[[ "$DATASET" == "all" ]] && aim_run 13_cross_dataset_keys
aim_run 14_normalization_diagnostics
aim_run 15_sample_relationships
aim_run 16_expression_patterns
aim_run 17_gsea_panels

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
