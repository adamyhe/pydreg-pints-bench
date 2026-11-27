#!/bin/bash
# Phase 2 : pydreg (GPU) + PINTS (CPU) on the 5 hg38 libraries, then phase-2 labels, then
# evaluation + bootstrap over ALL 10 libraries (phase-1 libraries pick up the new K562 PINTS-on-PRO-cap
# labels), then the report. Data must be fetched first: bash scripts/00b_fetch_phase2.sh
set -euo pipefail
p() { bash "$(dirname "${BASH_SOURCE[0]}")/sub.sh" "$@"; }
D=$(p --export=ALL,PPB_PHASE=2 slurm/10_pydreg.sbatch)
P=$(p --export=ALL,PPB_PHASE=2 slurm/11_pints.sbatch)
L=$(p --dependency=afterok:$P slurm/21_labels_phase2.sbatch)
E=$(p --export=ALL,PPB_PHASE=all --array=1-10 --dependency=afterok:$L:$D slurm/30_evaluate.sbatch)
B=$(p --export=ALL,PPB_PHASE=all --array=1-10 --dependency=afterok:$L:$D slurm/40_bootstrap.sbatch)
R=$(p --dependency=afterok:$E:$B slurm/41_report.sbatch)
echo "pydreg=$D pints=$P labels2=$L evaluate=$E bootstrap=$B report=$R"
