#!/bin/bash -l
# Evaluate each library (needs pydreg + PINTS outputs and built labels).
#SBATCH --partition=regular,long7,long30
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=1-5
#SBATCH --job-name=ppb_eval

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
LIB=$(lib_for_task "$SLURM_ARRAY_TASK_ID")
export OMP_NUM_THREADS="$SLURM_NTASKS"
cd "$PROJ"
"$VENV/python" scripts/03_evaluate.py "$LIB"
