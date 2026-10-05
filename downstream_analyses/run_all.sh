#!/usr/bin/env bash
# Run every step-3 analysis in order, then build REPORT.md.
#
# Reads csv_exports/ only, so it runs on a laptop in about a minute. Needs
# numpy, pandas, scipy (matplotlib optional: without it, figures are skipped).
#
#   ./downstream_analyses/run_all.sh
#   PYTHON=/path/to/python ./downstream_analyses/run_all.sh
#   AI_ADATA_OUT_ROOT=/elsewhere/csv_exports ./downstream_analyses/run_all.sh
#
# Order matters: 02 reads 01's sex calls; 04, 06 and 09 read 03's trends; 09 reads 08's modules.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"

for script in "$HERE"/[0-9][0-9]_*.py; do
    echo "=================================================================="
    echo ">>> $(basename "$script")"
    "$PYTHON" "$script"
done
echo "=================================================================="
"$PYTHON" "$HERE/build_report.py"
