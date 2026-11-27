#!/bin/bash
# Fetch everything the pipeline needs into data/ (and kentUtils into tools/kent if needed). Login node with
# internet; idempotent, so re-run it after an interruption. About 4 GB in total, most of it bigWigs and the
# two genome .2bit files. Logs go to logs/.
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$PROJ/logs"
bash "$PROJ/scripts/00_setup_tools.sh"
bash "$PROJ/scripts/00_fetch_data.sh"   2>&1 | tee "$PROJ/logs/00_fetch_data.log"
bash "$PROJ/scripts/00b_fetch_phase2.sh" 2>&1 | tee "$PROJ/logs/00b_fetch_phase2.log"
bash "$PROJ/scripts/00c_fetch_jurkat.sh" 2>&1 | tee "$PROJ/logs/00c_fetch_jurkat.log"
echo "all inputs in $PROJ/data"
