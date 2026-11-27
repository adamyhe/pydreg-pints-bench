#!/bin/bash -l
# Build the PINTS-paper and dREG-paper label sets (tool-independent; needs only downloaded data).
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=ppb_labels

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
cd "$PROJ"
"$VENV/python" scripts/01_build_labels.py
