#!/bin/bash
# Re-evaluate the 10 phase-1/2 libraries (adds dhs_tss_silent negatives) and re-bootstrap all 15 libraries
# (adds GC-matched and CpG-island-split conditions), then the report.
# Optional arg: colon-separated job IDs to wait for first (e.g. the phase-3 report job).
set -euo pipefail
cd /home2/ayh8/pydreg_pints
p() { sbatch --parsable "$@" | cut -d';' -f1; }
DEP=${1:+--dependency=afterany:$1}
E=$(p --export=ALL,PPB_PHASE=all --array=1-10 $DEP slurm/30_evaluate.sbatch)
B=$(p --export=ALL,PPB_PHASE=all --array=1-15 $DEP slurm/40_bootstrap.sbatch)
R=$(p --dependency=afterany:$E:$B slurm/41_report.sbatch)
echo "evaluate=$E bootstrap=$B report=$R"
