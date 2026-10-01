# Shared setup for the ai-adata-miner Slurm jobs. Sourced, never submitted.
#
# The #SBATCH header and `module load` stay inline in each job script (Slurm
# only reads directives from the submitted file, and keeping the module line
# visible there matches the rest of running_scripts/). What lives here is the
# part that would otherwise be copied five times: environment resolution, the
# per-step runner, and the end-of-job report.

AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
AIM_ENV="${AIM_ENV:-/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new}"
AIM_JOBS_OUT="${AIM_JOBS_OUT:-/miridan-data/annaludmir/jobs_output}"
MAILTO="${MAILTO:-annaludmir@mail.tau.ac.il}"

# Compute nodes accept mail (sendmail exits 0) but do not relay it, so the
# summary file is the reliable copy -- same finding as python_cell_cycle_annotation.sh.
# Flip to "true" only if a test message from a compute node actually arrives;
# otherwise use mail_job_summary.sh from the login node.
EMAIL_OUTPUT="${EMAIL_OUTPUT:-false}"

JOBID="${SLURM_JOB_ID:-local}"
OUTFILE="${AIM_JOBS_OUT}/${JOBID}.out"
SUMMARY="${AIM_JOBS_OUT}/${JOBID}_summary.txt"

# Pipeline configuration, read by config.py.
export AI_ADATA_DATA_ROOT="${AI_ADATA_DATA_ROOT:-/miridan-data/annaludmir/ndd_gene_modules/data}"
export AI_ADATA_OUT_ROOT="${AI_ADATA_OUT_ROOT:-${AIM_ROOT}/csv_exports}"

# Which dataset(s) to process, and the knobs worth changing per submission.
DATASET="${DATASET:-all}"
CHUNK_SIZE="${CHUNK_SIZE:-50000}"
TOP_GENES="${TOP_GENES:-12000}"
LIMIT_CELLS="${LIMIT_CELLS:-}"      # set to e.g. 20000 for a smoke test

AIM_FAILED=()
AIM_OK=()
AIM_START=$SECONDS


aim_setup() {
  cd "$AIM_ROOT"
  mkdir -p "$AIM_JOBS_OUT"
  echo "=============================================================="
  echo "job        : ${JOBID} on $(hostname)"
  echo "repo       : ${AIM_ROOT}"
  echo "env        : ${AIM_ENV}"
  echo "data root  : ${AI_ADATA_DATA_ROOT}"
  echo "csv exports: ${AI_ADATA_OUT_ROOT}"
  echo "dataset    : ${DATASET}"
  echo "chunk size : ${CHUNK_SIZE}"
  [[ -n "$LIMIT_CELLS" ]] && echo "LIMIT_CELLS: ${LIMIT_CELLS}  (SMOKE TEST -- results are partial)"
  echo "started    : $(date '+%F %T')"
  echo "=============================================================="
}


# aim_run <script_stem> [extra args...]
# Runs one pipeline script. Records failure and carries on, so one broken step
# does not throw away the steps that already succeeded -- important when step 09
# has already spent hours streaming the matrix.
aim_run() {
  local stem="$1"; shift
  local args=(--dataset "$DATASET" --chunk-size "$CHUNK_SIZE")
  [[ -n "$LIMIT_CELLS" ]] && args+=(--limit-cells "$LIMIT_CELLS")
  args+=("$@")

  echo ""
  echo "--------------------------------------------------------------"
  echo ">>> ${stem}   [$(date '+%T')]"
  echo "--------------------------------------------------------------"
  local t0=$SECONDS
  if mamba run -p "$AIM_ENV" python -u "scripts_generated/${stem}.py" "${args[@]}"; then
    local dt=$(( SECONDS - t0 ))
    echo "--- ${stem} OK (${dt}s)"
    AIM_OK+=("${stem} ${dt}s")
  else
    local rc=$? dt=$(( SECONDS - t0 ))
    echo "!!! ${stem} FAILED rc=${rc} after ${dt}s -- continuing"
    AIM_FAILED+=("${stem} rc=${rc}")
  fi
}


aim_csv_count() {
  find "$AI_ADATA_OUT_ROOT" -name '*.csv' 2>/dev/null | wc -l | tr -d ' '
}


# End-of-job report. Slurm's own --mail-type mail carries only the job status,
# never the output, so the interesting part of the log is distilled here.
aim_report() {
  rc=$?                          # must stay first: anything else overwrites it
  local elapsed=$(( SECONDS - AIM_START ))
  {
    echo "job ${JOBID}  host=$(hostname)  exit=${rc}  elapsed=$((elapsed/3600))h$(( (elapsed%3600)/60 ))m"
    echo "stage: ${AIM_STAGE:-unnamed}   dataset: ${DATASET}"
    echo "log:   ${OUTFILE}"
    [[ -n "$LIMIT_CELLS" ]] && echo "NOTE:  LIMIT_CELLS=${LIMIT_CELLS} -- this was a smoke test, not a full run"
    echo
    echo "steps completed (${#AIM_OK[@]}):"
    printf '  %s\n' "${AIM_OK[@]:-none}"
    if (( ${#AIM_FAILED[@]} )); then
      echo
      echo "STEPS FAILED (${#AIM_FAILED[@]}):"
      printf '  %s\n' "${AIM_FAILED[@]}"
    fi
    echo
    echo "CSV files now under ${AI_ADATA_OUT_ROOT}: $(aim_csv_count)"
    for m in "${AI_ADATA_OUT_ROOT}"/*/_manifest.csv; do
      [[ -e "$m" ]] && echo "  manifest $(basename "$(dirname "$m")"): $(( $(wc -l < "$m") - 1 )) entries"
    done

    # Surface the findings that change how the results must be read, rather
    # than making someone grep a multi-thousand-line log for them.
    if [[ -f "$OUTFILE" ]]; then
      echo
      echo "--- warnings and flagged findings in the log ---"
      # Exclude the "wrote <file>  <description>" echoes first: those carry the
      # same words in their descriptions and would bury the real findings.
      grep -v '  wrote ' "$OUTFILE" \
        | grep -iE 'WARNING|ERROR|strongest entanglement|dims flagged as likely technical|'\
'exceeds|skipping|no varm|no spliced|enrichments at FDR|widest TMM spread|'\
'groups show a global shift|K maximising|below .* cells|not found' \
        | sed 's/^\[[0-9:]*\] *//' | sort -u | head -40 || echo "(none)"
    fi
  } > "$SUMMARY" 2>/dev/null || true

  if [[ "$EMAIL_OUTPUT" == "true" && -x /usr/sbin/sendmail ]]; then
    {
      echo "To: ${MAILTO}"
      echo "Subject: [ai-adata-miner/${AIM_STAGE:-run}] job ${JOBID} exit=${rc}"
      echo "From: ${MAILTO}"
      echo
      cat "$SUMMARY"
    } | /usr/sbin/sendmail -t || true
  fi

  echo ""
  echo "=============================================================="
  cat "$SUMMARY" 2>/dev/null || true
}
