#!/bin/bash -l
# Paired bootstrap CIs per library (recomputes element scores; needs eval inputs).
#SBATCH --clusters=cbsubscb
#SBATCH --account=danko_0001
#SBATCH --partition=regular,long7,long30
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-5
#SBATCH --job-name=ppb_boot
#SBATCH --output=/home/ayh8/slurm/ppb_boot.%A_%a.out

source /home2/ayh8/pydreg_pints/slurm/common.sh
make_scratch
LIB=$(lib_for_task "$SLURM_ARRAY_TASK_ID")
export OMP_NUM_THREADS="$SLURM_NTASKS"
cd "$PROJ"
"$VENV/python" scripts/05_bootstrap.py --n-boot 1000 "$LIB"
