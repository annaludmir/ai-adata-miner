#!/usr/bin/env bash
# Run every step-3 analysis in order, then build REPORT.md.
#
# Reads csv_exports/ only, so it runs on a laptop in a few minutes. Needs
# numpy, pandas, scipy (matplotlib optional: without it, figures are skipped).
#
#   ./downstream_analyses/run_all.sh
#   PYTHON=/path/to/python ./downstream_analyses/run_all.sh
#   AI_ADATA_OUT_ROOT=/elsewhere/csv_exports ./downstream_analyses/run_all.sh
#
# Order matters: 02 and 14 read 01's sex calls; 04, 06, 09, 10 and 19 read 03's trends; 09, 10,
# 13, 15, 18 and 19 use 08's modules. Number order satisfies all of these.
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
