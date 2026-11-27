#!/bin/bash -l
# Phase-3 (Jurkat) labels (needs PINTS calls on Jurkat_CAGE).
#SBATCH --clusters=cbsubscb
#SBATCH --account=danko_0001
#SBATCH --partition=short
#SBATCH --ntasks=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --job-name=ppb_labels3
#SBATCH --output=/home/ayh8/slurm/ppb_labels3.%j.out

source /home2/ayh8/pydreg_pints/slurm/common.sh
make_scratch
cd "$PROJ"
"$VENV/python" scripts/01c_build_labels_jurkat.py
