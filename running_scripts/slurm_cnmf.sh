#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_cnmf
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# cNMF: gene programmes with no prior (script 26), one stratum per job -- the
# prior-free counterpart of the Spectra fits. Uses the Spectra subsample
# ($AIM_WORK/spectra/<ns>/input.h5ad; script 24 runs first if it is missing)
# and writes csv_exports/<ns>/26_cnmf/. Step-3 analysis 36 reads it.
#
# 7 values of K x 100 NMF replicates, spread over the job's CPUs (~2-4 h).
# The env needs:  /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new/bin/python -m pip install cnmf
#
#   for ds in cortex human_dev; do for ch in v2 v3; do
#     sbatch --export=ALL,DATASET=$ds,CHEMISTRY=$ch slurm_cnmf.sh; done; done
# Another K from the saved replicates (consensus only, minutes): CNMF_K=30 CNMF_SKIP_FACTORIZE=true

set -uo pipefail

AIM_STAGE="cnmf"
AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
source "${AIM_ROOT}/running_scripts/_common.sh"
trap aim_report EXIT

export AIM_WORK="${AIM_WORK:-/miridan-data/annaludmir/aim_work}"
CNMF_SKIP_FACTORIZE="${CNMF_SKIP_FACTORIZE:-false}"
[[ "$DATASET" == "all" ]] && { echo "set DATASET=cortex or human_dev (one stratum per job)"; exit 1; }
[[ "$CHEMISTRY" == "each" ]] && { echo "set CHEMISTRY=v2 or v3 (one stratum per job)"; exit 1; }

module load mamba/mamba-1.5.8
mamba activate "$AIM_ENV"
aim_setup
echo "work dir   : ${AIM_WORK}"
echo "K          : ${CNMF_K:-chosen by stability}"
if ! mamba run -p "$AIM_ENV" python -c "import cnmf" 2>/dev/null; then
  echo "ERROR: cnmf not importable in $AIM_ENV -- install it:"
  echo "  $AIM_ENV/bin/python -m pip install cnmf"
  exit 1
fi

[[ -f "${AIM_WORK}/spectra/${DATASET}__${CHEMISTRY}/input.h5ad" ]] || aim_run 24_spectra_input
if [[ "$CNMF_SKIP_FACTORIZE" == "true" ]]; then
  aim_run 26_cnmf --skip-factorize
else
  aim_run 26_cnmf
fi

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
