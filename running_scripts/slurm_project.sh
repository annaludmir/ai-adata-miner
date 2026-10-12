#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_project
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Score every cell with the fitted gene programmes (script 28), one stratum per
# job. Holds each fit's gene side fixed (Spectra full / core / STRING, cNMF,
# scHPF -- whichever exist for the stratum) and estimates every cell's
# activities from the stratum's h5ad (streamed), so analysis 36's age, region
# and donor trends rest on all cells. Writes csv_exports/<ns>/28_projection/<fit>/
# incl. projection_check.csv (projected vs fitted on the subsample: r near 1).
#
# Run after the fits; ~1-2 h per human_dev stratum (Spectra's Poisson fit dominates).
#   for ds in cortex human_dev; do for ch in v2 v3; do
#     sbatch --export=ALL,DATASET=$ds,CHEMISTRY=$ch slurm_project.sh; done; done
# Only some fits: PROJECT_FITS="cnmf schpf"

set -uo pipefail

AIM_STAGE="project"
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

export AIM_WORK="${AIM_WORK:-/miridan-data/annaludmir/aim_work}"
[[ "$DATASET" == "all" ]] && { echo "set DATASET=cortex or human_dev (one stratum per job)"; exit 1; }
[[ "$CHEMISTRY" == "each" ]] && { echo "set CHEMISTRY=v2 or v3 (one stratum per job)"; exit 1; }

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup
echo "work dir   : ${AIM_WORK}"
echo "fits       : ${PROJECT_FITS:-all that exist}"

if [[ -n "${PROJECT_FITS:-}" ]]; then
  # shellcheck disable=SC2086
  aim_run 28_project_programs --fits $PROJECT_FITS
else
  aim_run 28_project_programs
fi

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
