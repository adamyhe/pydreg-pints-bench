#!/bin/bash -l
# Phase-2 labels (needs PINTS calls on the PRO-cap libraries).
#SBATCH --clusters=cbsubscb
#SBATCH --account=danko_0001
#SBATCH --partition=short
#SBATCH --ntasks=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --job-name=ppb_labels2
#SBATCH --output=/home/ayh8/slurm/ppb_labels2.%j.out

source /home2/ayh8/pydreg_pints/slurm/common.sh
make_scratch
cd "$PROJ"
"$VENV/python" scripts/01b_build_labels_phase2.py
