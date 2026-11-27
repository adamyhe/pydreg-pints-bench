#!/bin/bash
# Memory profiles of PINTS (6 tasks) and pydreg (2 tasks); results in results/memprof/*.mem.{json,tsv}.
set -euo pipefail
p() { bash "$(dirname "${BASH_SOURCE[0]}")/sub.sh" "$@"; }
P=$(p slurm/50_memprof_pints.sbatch)
D=$(p slurm/51_memprof_pydreg.sbatch)
echo "memprof pints=$P pydreg=$D  ->  results/memprof/"
