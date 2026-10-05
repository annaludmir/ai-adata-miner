#!/usr/bin/env bash
# Run the full extraction pipeline (steps 1 & 2) for both datasets.
#
#   ./run_all.sh                     # everything, both datasets
#   ./run_all.sh --dataset cortex    # one dataset
#   ./run_all.sh --limit-cells 20000 # smoke test on a slice first (recommended)
#   ./run_all.sh --chemistry all     # pool v2+v3 (default is to run each apart)
#   TOP_GENES=20000 ./run_all.sh     # env var: only step 09 accepts --top-genes
#
# Scripts 00-08 and 13 read only .obs/.var/.obsm/.varm and finish in minutes.
# Script 09 streams X and is the long one; 10-11 then work purely on its CSVs.
set -euo pipefail

cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
ARGS=("$@")

ORDER=(
  00_inventory_and_schema      # ground-truth inventory; run first
  01_obs_overview_and_qc
  02_composition
  03_cellcycle
  04_cluster_profiles
  05_covariate_confounds
  06_embeddings
  07_factor_modules            # needs varm/Loadings (human_dev only)
  08_factor_activity
  09_pseudobulk                # streams X -- the expensive step
  10_marker_specificity        # reads 09's CSVs
  11_gene_panels               # reads 09's CSVs
  12_splicing_layers           # needs spliced/unspliced layers (cortex only)
  13_cross_dataset_keys        # needs both datasets
  14_normalization_diagnostics # TMM factors + MA diagnostics, reads 09's CSVs
  15_sample_relationships      # correlation / hierarchical clustering / PCA
  16_expression_patterns       # Z-scored K-means expression patterns
  17_gsea_panels               # GSEA of the panels, reads 10's CSVs
  18_chemistry_comparability   # pooled by design: what survives stratification
  19_rg_subtypes               # oRG vs vRG per radial-glia cell (streams X; after 09)
)

mkdir -p logs
FAILED=()
for script in "${ORDER[@]}"; do
  echo ""
  echo "=============================================================="
  echo ">>> ${script}"
  echo "=============================================================="
  # --top-genes is only accepted by 09; passing it to the others is an error,
  # so it is added per script rather than forwarded blindly.
  STEP_ARGS=("${ARGS[@]}")
  if [[ "$script" == 09_* && -n "${TOP_GENES:-}" ]]; then
    STEP_ARGS+=(--top-genes "$TOP_GENES")
  fi
  if "$PYTHON" "scripts_generated/${script}.py" "${STEP_ARGS[@]}" 2>&1 | tee "logs/${script}.log"; then
    echo "--- ${script} OK"
  else
    echo "!!! ${script} FAILED (see logs/${script}.log) -- continuing"
    FAILED+=("${script}")
  fi
done

echo ""
echo "=============================================================="
# Honour AI_ADATA_OUT_ROOT rather than assuming the repo-local folder.
OUT_ROOT="$("$PYTHON" -c 'import config; print(config.CSV_EXPORTS)')"
echo "CSV output root: ${OUT_ROOT}"
echo "CSV files written: $(find "$OUT_ROOT" -name '*.csv' 2>/dev/null | wc -l | tr -d ' ')"
for m in "$OUT_ROOT"/*/_manifest.csv; do
  [ -e "$m" ] && echo "  manifest: $m ($(( $(wc -l < "$m") - 1 )) entries)"
done
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "FAILED: ${FAILED[*]}"
  exit 1
fi
echo "All steps completed."
