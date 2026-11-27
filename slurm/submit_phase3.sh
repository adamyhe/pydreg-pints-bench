#!/bin/bash
# Phase 3 (Jurkat, hg19) : pydreg (GPU) + PINTS on the 5 Jurkat libraries, Jurkat labels,
# evaluation + bootstrap of those libraries, then the report over everything.
# Data first: bash scripts/00c_fetch_jurkat.sh
# Optional arg: colon-separated job IDs the final report should also wait for (e.g. phase-2 eval:boot).
set -euo pipefail
p() { bash "$(dirname "${BASH_SOURCE[0]}")/sub.sh" "$@"; }
EXTRA=${1:+:$1}
D=$(p --export=ALL,PPB_PHASE=3 slurm/10_pydreg.sbatch)
P=$(p --export=ALL,PPB_PHASE=3 slurm/11_pints.sbatch)
L=$(p --dependency=afterok:$P slurm/22_labels_jurkat.sbatch)
E=$(p --export=ALL,PPB_PHASE=3 --array=1-5 --dependency=afterok:$L:$D slurm/30_evaluate.sbatch)
B=$(p --export=ALL,PPB_PHASE=3 --array=1-5 --dependency=afterok:$L:$D slurm/40_bootstrap.sbatch)
R=$(p --dependency=afterany:$E:$B$EXTRA slurm/41_report.sbatch)
echo "pydreg=$D pints=$P labels3=$L evaluate=$E bootstrap=$B report=$R"
