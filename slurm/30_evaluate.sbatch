#!/bin/bash -l
# Evaluate each library (needs pydreg + PINTS outputs and built labels).
#SBATCH --clusters=cbsubscb
#SBATCH --account=danko_0001
#SBATCH --partition=regular,long7,long30
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=1-5
#SBATCH --job-name=ppb_eval
#SBATCH --output=/home/ayh8/slurm/ppb_eval.%A_%a.out

source /home2/ayh8/pydreg_pints/slurm/common.sh
make_scratch
LIB=$(lib_for_task "$SLURM_ARRAY_TASK_ID")
export OMP_NUM_THREADS="$SLURM_NTASKS"
cd "$PROJ"
"$VENV/python" scripts/03_evaluate.py "$LIB"
