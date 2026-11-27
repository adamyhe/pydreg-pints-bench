#!/bin/bash -l
# Phase-2 labels (needs PINTS calls on the PRO-cap libraries).
#SBATCH --partition=short
#SBATCH --ntasks=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --job-name=ppb_labels2

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
cd "$PROJ"
"$VENV/python" scripts/01b_build_labels_phase2.py
