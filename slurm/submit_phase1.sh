#!/bin/bash
# Submit phase 1  with dependencies. Run from a login node AFTER:
#   1) bash scripts/00_fetch_data.sh   (downloads)
#   2) sbatch slurm/00_gpu_smoke.sbatch and check its log
set -euo pipefail
p() { bash "$(dirname "${BASH_SOURCE[0]}")/sub.sh" "$@"; }
L=$(p slurm/20_labels.sbatch)
D=$(p slurm/10_pydreg.sbatch)
P=$(p slurm/11_pints.sbatch)
E=$(p --dependency=afterok:$L:$D:$P slurm/30_evaluate.sbatch)
echo "labels=$L pydreg=$D pints=$P evaluate=$E"
echo "after evaluate finishes: .venv/bin/python scripts/04_summarize.py"
