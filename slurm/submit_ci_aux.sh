#!/bin/bash
# Re-evaluate (adds the dreg_aux label family) and bootstrap CIs in parallel, then the report.
# Labels must already exist under data/labels/built (dreg_aux built 2026-11-26).
set -euo pipefail
p() { bash "$(dirname "${BASH_SOURCE[0]}")/sub.sh" "$@"; }
E=$(p slurm/30_evaluate.sbatch)
B=$(p slurm/40_bootstrap.sbatch)
R=$(p --dependency=afterok:$E:$B slurm/41_report.sbatch)
echo "evaluate=$E bootstrap=$B report=$R  ->  results/summary/report.md"
