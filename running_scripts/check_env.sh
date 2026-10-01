#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=FAIL
#SBATCH --job-name=aim_checkenv
#SBATCH --mem=32G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:20:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Does AIM_ENV actually satisfy this pipeline? Run this before anything else.
#
# Light enough for the login node, and submittable if you would rather not run
# anything there:
#
#   bash check_env.sh            # login node, ~1 minute
#   sbatch check_env.sh          # as a job
#   AIM_ENV=/path/to/other/env bash check_env.sh
#
# It checks the five packages the pipeline needs, resolves anndata's read_elem
# (which has moved twice between versions), imports the repo's own modules, and
# then does the real test: opens both h5ad files and reads their shape and a
# slice of .obs. An env that passes this can run the pipeline.

set -uo pipefail

AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
AIM_ENV="${AIM_ENV:-/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new}"
export AI_ADATA_DATA_ROOT="${AI_ADATA_DATA_ROOT:-/miridan-data/annaludmir/ndd_gene_modules/data}"

module load mamba/mamba-1.5.8 2>/dev/null || true
mamba activate "$AIM_ENV" 2>/dev/null || true
cd "$AIM_ROOT" || { echo "FAIL: AIM_ROOT=$AIM_ROOT does not exist"; exit 1; }

echo "env : $AIM_ENV"
echo "repo: $AIM_ROOT"
echo "data: $AI_ADATA_DATA_ROOT"
echo

mamba run -p "$AIM_ENV" python -u "${AIM_ROOT}/running_scripts/check_env.py"
rc=$?
echo "Python exit code: $rc"
exit $rc
