#!/bin/bash
# submit_all.sh -- submit the pipeline as a dependency chain. RUN ON THE LOGIN NODE.
#
#   ./submit_all.sh                      # full run, both datasets
#   LIMIT_CELLS=20000 ./submit_all.sh    # smoke test first (recommended)
#   DATASETS="cortex" ./submit_all.sh    # one dataset only
#   MAIL_SUMMARY=true ./submit_all.sh    # mail the final summary when it lands
#   DOWNSTREAM=false ./submit_all.sh     # stop after stage 3 (skip step-3 analyses)
#
# Shape of the chain:
#
#   stage 1  metadata (both datasets, 16G)
#      |
#      +--> stage 2  pseudobulk cortex    (128G)  ─┐  run concurrently: they
#      +--> stage 2  pseudobulk human_dev (128G)  ─┤  read different files and
#      |                                           │  write different folders
#      v (afterok, both)
#   stage 3  derived analyses (16G, CSV only)
#      |
#      v (afterok)
#   stage 4  downstream_analyses/ + REPORT.md (8G), written outside the repo
#            to AIM_DOWNSTREAM_OUT (default /miridan-data/annaludmir/aim_downstream)
#
# Splitting stage 2 per dataset is the point of doing this rather than
# slurm_full_pipeline.sh: the two matrix passes are independent, and human_dev
# is ~5x the cells of cortex, so running them together roughly halves wall time.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

DATASETS="${DATASETS:-cortex human_dev}"
# 'each' runs v2 and v3 separately. Stage 2 fans out over dataset x chemistry so
# the extra pass costs cluster slots rather than wall time.
CHEMISTRY="${CHEMISTRY:-each}"
MAIL_SUMMARY="${MAIL_SUMMARY:-false}"
DOWNSTREAM="${DOWNSTREAM:-true}"
MAIL_HELPER="${MAIL_HELPER:-/miridan-data/annaludmir/ndd_gene_modules/running_scripts/mail_job_summary.sh}"

# Everything the job scripts read from the environment has to be forwarded
# explicitly; --export=ALL alone does not carry through a dependency chain.
pass_env() {
  local kv="ALL"
  for v in AIM_ROOT AIM_ENV AI_ADATA_DATA_ROOT AI_ADATA_OUT_ROOT \
           CHUNK_SIZE TOP_GENES LIMIT_CELLS AIM_EXCLUSIONS AIM_DOWNSTREAM_OUT \
           EMAIL_OUTPUT; do
    [[ -n "${!v:-}" ]] && kv="${kv},${v}=${!v}"
  done
  echo "$kv"
}

[[ -n "${LIMIT_CELLS:-}" ]] && \
  echo "NOTE: LIMIT_CELLS=${LIMIT_CELLS} -- this is a SMOKE TEST; outputs are partial." && echo

STAGE1=$(sbatch --parsable --export="$(pass_env),DATASET=all,CHEMISTRY=${CHEMISTRY}" slurm_01_metadata.sh)
printf 'stage 1  %-24s -> job %s\n' "metadata (all)" "$STAGE1"

case "$CHEMISTRY" in
  each) CHEMS="v2 v3" ;;
  all)  CHEMS="all" ;;
  *)    CHEMS="$CHEMISTRY" ;;
esac

STAGE2_IDS=()
for ds in $DATASETS; do
  for ch in $CHEMS; do
    jid=$(sbatch --parsable --dependency=afterok:"${STAGE1}" \
                 --job-name="aim_pbulk_${ds}_${ch}" \
                 --export="$(pass_env),DATASET=${ds},CHEMISTRY=${ch}" \
                 slurm_02_pseudobulk.sh)
    STAGE2_IDS+=("$jid")
    printf 'stage 2  %-24s -> job %s  (after %s)\n' "pseudobulk ${ds} ${ch}" "$jid" "$STAGE1"
  done
done

DEP=$(IFS=:; echo "${STAGE2_IDS[*]}")
STAGE3=$(sbatch --parsable --dependency=afterok:"${DEP}" \
                --export="$(pass_env),DATASET=all,CHEMISTRY=${CHEMISTRY}" slurm_03_derived.sh)
printf 'stage 3  %-24s -> job %s  (after %s)\n' "derived analyses (all)" "$STAGE3" "$DEP"
LAST="$STAGE3"
STAGE4=""
if [[ "$DOWNSTREAM" == "true" ]]; then
  STAGE4=$(sbatch --parsable --dependency=afterok:"${STAGE3}" \
                  --export="$(pass_env),DATASET=all,CHEMISTRY=${CHEMISTRY}" slurm_04_downstream.sh)
  printf 'stage 4  %-24s -> job %s  (after %s)\n' "downstream analyses" "$STAGE4" "$STAGE3"
  LAST="$STAGE4"
fi

echo
echo "watch:    squeue -u \$USER"
ALL_IDS=("$STAGE1" "${STAGE2_IDS[@]}" "$STAGE3" ${STAGE4:+"$STAGE4"})
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
    nohup "$MAIL_HELPER" "$LAST" >/dev/null 2>&1 &
    echo
    echo "mail_job_summary.sh watching job ${LAST} (compute nodes cannot relay mail themselves)"
  else
    echo
    echo "WARNING: MAIL_SUMMARY=true but ${MAIL_HELPER} is not executable -- skipping."
  fi
fi
