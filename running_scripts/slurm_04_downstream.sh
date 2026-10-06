#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_downstream
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Stage 4 -- step-3 downstream analyses (downstream_analyses/), CSV-only.
#
# Results go OUTSIDE the repo by default (AIM_DOWNSTREAM_OUT): REPORT.md, the
# SUMMARY.md files and the figures are tracked in git, so writing them into
# the cluster checkout would make the next `git pull` conflict. Copy the folder
# back with rsync, or point AIM_DOWNSTREAM_OUT elsewhere.
#
# Analyses 01-05 take ~20 s on a laptop; the co-expression analyses 06-08 add
# about a minute and ~1.5 GB (longer with large gene lists). Cluster nodes ran
# stage 3 several times slower than a laptop, hence 8G / 1 CPU / 1h.
#
#   sbatch slurm_04_downstream.sh
#   sbatch --dependency=afterok:<stage3_jobid> slurm_04_downstream.sh
#   AIM_DOWNSTREAM_OUT=/elsewhere/results sbatch --export=ALL slurm_04_downstream.sh

set -uo pipefail

AIM_STAGE="04_downstream"
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

export AIM_DOWNSTREAM_OUT="${AIM_DOWNSTREAM_OUT:-/miridan-data/annaludmir/aim_downstream/results}"
# Analysis 28 compares these results with a strict-QC run (AIM_STRICT=true) when one exists.
export AIM_STRICT_RESULTS="${AIM_STRICT_RESULTS:-/miridan-data/annaludmir/aim_downstream_strict/results}"
mkdir -p "$AIM_DOWNSTREAM_OUT"

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup
echo "step-3 out : ${AIM_DOWNSTREAM_OUT}  (REPORT.md one level up)"

if ! ls "${AI_ADATA_OUT_ROOT}"/*/09_pseudobulk/*__pseudobulk_counts.csv >/dev/null 2>&1; then
  echo "ERROR: no pseudobulk CSVs under ${AI_ADATA_OUT_ROOT}."
  echo "       Run stages 1-3 (submit_all.sh) first."
  exit 1
fi

# In number order: 02 reads 01's sex calls, 04 reads 03's trends.
for script in downstream_analyses/[0-9][0-9]_*.py; do
  aim_run_script "$(basename "$script" .py)" "$script"
done
aim_run_script build_report downstream_analyses/build_report.py

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
