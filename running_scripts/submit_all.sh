#!/bin/bash
# submit_all.sh -- submit the pipeline as a dependency chain. RUN ON THE LOGIN NODE.
#
#   ./submit_all.sh                      # full run, both datasets
#   LIMIT_CELLS=20000 ./submit_all.sh    # smoke test first (recommended)
#   DATASETS="cortex" ./submit_all.sh    # one dataset only
#   MAIL_SUMMARY=true ./submit_all.sh    # mail the final summary when it lands
#
# Shape of the chain:
#
#   stage 1  metadata (both datasets, 250G)
#      |
#      +--> stage 2  pseudobulk cortex    (500G)  ─┐  run concurrently: they
#      +--> stage 2  pseudobulk human_dev (500G)  ─┤  read different files and
#      |                                           │  write different folders
#      v (afterok, both)
#   stage 3  derived analyses (200G, CSV only)
#
# Splitting stage 2 per dataset is the point of doing this rather than
# slurm_full_pipeline.sh: the two matrix passes are independent, and human_dev
# is ~5x the cells of cortex, so running them together roughly halves wall time.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

DATASETS="${DATASETS:-cortex human_dev}"
MAIL_SUMMARY="${MAIL_SUMMARY:-false}"
MAIL_HELPER="${MAIL_HELPER:-/miridan-data/annaludmir/ndd_gene_modules/running_scripts/mail_job_summary.sh}"

# Everything the job scripts read from the environment has to be forwarded
# explicitly; --export=ALL alone does not carry through a dependency chain.
pass_env() {
  local kv="ALL"
  for v in AIM_ROOT AIM_ENV AI_ADATA_DATA_ROOT AI_ADATA_OUT_ROOT \
           CHUNK_SIZE TOP_GENES LIMIT_CELLS EMAIL_OUTPUT; do
    [[ -n "${!v:-}" ]] && kv="${kv},${v}=${!v}"
  done
  echo "$kv"
}

[[ -n "${LIMIT_CELLS:-}" ]] && \
  echo "NOTE: LIMIT_CELLS=${LIMIT_CELLS} -- this is a SMOKE TEST; outputs are partial." && echo

STAGE1=$(sbatch --parsable --export="$(pass_env),DATASET=all" slurm_01_metadata.sh)
printf 'stage 1  %-24s -> job %s\n' "metadata (all)" "$STAGE1"

STAGE2_IDS=()
for ds in $DATASETS; do
  jid=$(sbatch --parsable --dependency=afterok:"${STAGE1}" \
               --job-name="aim_pbulk_${ds}" \
               --export="$(pass_env),DATASET=${ds}" slurm_02_pseudobulk.sh)
  STAGE2_IDS+=("$jid")
  printf 'stage 2  %-24s -> job %s  (after %s)\n' "pseudobulk ${ds}" "$jid" "$STAGE1"
done

DEP=$(IFS=:; echo "${STAGE2_IDS[*]}")
STAGE3=$(sbatch --parsable --dependency=afterok:"${DEP}" \
                --export="$(pass_env),DATASET=all" slurm_03_derived.sh)
printf 'stage 3  %-24s -> job %s  (after %s)\n' "derived analyses (all)" "$STAGE3" "$DEP"

echo
echo "watch:    squeue -u \$USER"
ALL_IDS=("$STAGE1" "${STAGE2_IDS[@]}" "$STAGE3")
echo "logs:     ${AIM_JOBS_OUT:-/miridan-data/annaludmir/jobs_output}/{$(IFS=,; echo "${ALL_IDS[*]}")}.out"
echo "summary:  ${AIM_JOBS_OUT:-/miridan-data/annaludmir/jobs_output}/<jobid>_summary.txt"
echo "cancel:   scancel ${ALL_IDS[*]}"

# A failed stage leaves its dependents in the queue as DependencyNeverSatisfied;
# they never run, but Slurm does not always clear them. Worth knowing.
echo
echo "if a stage fails, its dependents stay queued as DependencyNeverSatisfied"
echo "-- clear them with the scancel line above."

if [[ "$MAIL_SUMMARY" == "true" ]]; then
  if [[ -x "$MAIL_HELPER" ]]; then
    nohup "$MAIL_HELPER" "$STAGE3" >/dev/null 2>&1 &
    echo
    echo "mail_job_summary.sh watching job ${STAGE3} (compute nodes cannot relay mail themselves)"
  else
    echo
    echo "WARNING: MAIL_SUMMARY=true but ${MAIL_HELPER} is not executable -- skipping."
  fi
fi
