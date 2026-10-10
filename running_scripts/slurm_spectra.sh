#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_spectra
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Spectra: knowledge-guided gene programmes (scripts 24 and 25), one stratum per job.
#
# 24 streams X once to draw a stratified ~25k-cell subsample into
#    $AIM_WORK/spectra/<dataset>__<chem>/input.h5ad
# 25 fits Spectra on it (prior: gene lists, seed panels, core GO processes,
#    class markers) and writes csv_exports/<ns>/25_spectra/.
# Needs stages 1-2 (and the annotation folder for the GO prior). Step-3
# analysis 36 reads the outputs.
#
# The conda env needs:  /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new/bin/python -m pip install scSpectra "setuptools<81"
# (torch comes with scSpectra; for the gpu backend install a CUDA build of torch).
#
# CPU (default; 2-6 h per stratum at 16 cores, up to 10000 epochs at ~2-4 s each):
#   for ds in cortex human_dev; do for ch in v2 v3; do
#     sbatch --export=ALL,DATASET=$ds,CHEMISTRY=$ch slurm_spectra.sh; done; done
# GPU (the package's minibatched Spectra_gpu module):
#   sbatch --partition=<gpu partition> --gres=gpu:1 \
#          --export=ALL,DATASET=human_dev,CHEMISTRY=v2,SPECTRA_BACKEND=gpu slurm_spectra.sh
# Refit without redrawing the subsample: SKIP_INPUT=true
# Core-prior fit (no gene lists / NDD seeds in the prior; -> 25_spectra_core/): SPECTRA_PRIOR=core

set -uo pipefail

AIM_STAGE="spectra"
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

export AIM_WORK="${AIM_WORK:-/miridan-data/annaludmir/aim_work}"
export SPECTRA_BACKEND="${SPECTRA_BACKEND:-cpu}"
SKIP_INPUT="${SKIP_INPUT:-false}"
[[ "$DATASET" == "all" ]] && { echo "set DATASET=cortex or human_dev (one stratum per job)"; exit 1; }
[[ "$CHEMISTRY" == "each" ]] && { echo "set CHEMISTRY=v2 or v3 (one stratum per job)"; exit 1; }

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup
echo "work dir   : ${AIM_WORK}"
echo "backend    : ${SPECTRA_BACKEND}${SPECTRA_EPOCHS:+, ${SPECTRA_EPOCHS} epochs}"
echo "prior      : ${SPECTRA_PRIOR:-full}"
if ! mamba run -p "$AIM_ENV" python -c "import Spectra" 2>/dev/null; then
  echo "ERROR: scSpectra not importable in $AIM_ENV -- install it:"
  echo "  mamba run -p \"$AIM_ENV\" pip install scSpectra \"setuptools<81\""
  exit 1
fi
[[ "$SPECTRA_BACKEND" == "gpu" ]] && mamba run -p "$AIM_ENV" python -c \
  "import torch; print('cuda available:', torch.cuda.is_available())"

[[ "$SKIP_INPUT" == "true" ]] || aim_run 24_spectra_input
aim_run 25_spectra_fit --backend "$SPECTRA_BACKEND"

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
