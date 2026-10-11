#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aim_schpf
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# scHPF: gene programmes with no prior (script 27), one stratum per job -- a
# second prior-free method beside cNMF. Uses the Spectra subsample
# ($AIM_WORK/spectra/<ns>/input.h5ad; script 24 runs first if it is missing)
# and K from the stratum's cNMF fit (run slurm_cnmf.sh, incl. its finalise step,
# first), and writes csv_exports/<ns>/27_schpf/. Step-3 analysis 36 reads it.
#
# The env needs scHPF from GitHub (it is not on PyPI):
#   /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new/bin/python -m pip install git+https://github.com/simslab/scHPF.git
#
#   for ds in cortex human_dev; do for ch in v2 v3; do
#     sbatch --export=ALL,DATASET=$ds,CHEMISTRY=$ch slurm_schpf.sh; done; done
# Other settings: SCHPF_K=30 (instead of cNMF's K), SCHPF_TRIALS=5

set -uo pipefail

AIM_STAGE="schpf"
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
echo "K          : ${SCHPF_K:-from the cNMF fit}"
if ! mamba run -p "$AIM_ENV" python -c "import schpf" 2>/dev/null; then
  echo "ERROR: scHPF not importable in $AIM_ENV -- install it:"
  echo "  $AIM_ENV/bin/python -m pip install git+https://github.com/simslab/scHPF.git"
  exit 1
fi

[[ -f "${AIM_WORK}/spectra/${DATASET}__${CHEMISTRY}/input.h5ad" ]] || aim_run 24_spectra_input
aim_run 27_schpf

(( ${#AIM_FAILED[@]} == 0 )) || exit 1
