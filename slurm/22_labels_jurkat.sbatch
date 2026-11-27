#!/bin/bash -l
# Phase-3 (Jurkat) labels (needs PINTS calls on Jurkat_CAGE).
#SBATCH --partition=short
#SBATCH --ntasks=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --job-name=ppb_labels3

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
cd "$PROJ"
"$VENV/python" scripts/01c_build_labels_jurkat.py
