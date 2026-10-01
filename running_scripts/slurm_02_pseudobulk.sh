#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_pbulk
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Stage 2 of 3 -- the only stage that reads the count matrix.
#
# Script 09 streams X once and fills every grouping's accumulator in that single
# pass; script 12 does the same over the spliced/unspliced layers (cortex only,
# human_dev has no layers and skips cleanly).
#
# Memory: a 50k-cell chunk of human_dev is roughly 2-3 GB sparse, and the pass
# holds about three copies of it (raw, CP10K-normalised, binarised) plus the
# dense per-group accumulators (~0.5 GB at 12k genes). 500G is generous; drop
# CHUNK_SIZE if a node ever comes up short.
#
#   sbatch slurm_02_pseudobulk.sh
#   DATASET=human_dev sbatch --export=ALL slurm_02_pseudobulk.sh
#   CHUNK_SIZE=25000 TOP_GENES=20000 sbatch --export=ALL slurm_02_pseudobulk.sh

set -uo pipefail

AIM_STAGE="02_pseudobulk"
# Slurm runs this from a spool copy on the compute node, so BASH_SOURCE does not
# point at the repo. Resolve the helper from an absolute root instead.
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup
echo "top genes  : ${TOP_GENES}"

aim_run 09_pseudobulk --top-genes "$TOP_GENES"
aim_run 12_splicing_layers     # needs spliced/unspliced layers -- cortex only

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
