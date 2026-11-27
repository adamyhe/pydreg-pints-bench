#!/bin/bash -l
# Build the PINTS-paper and dREG-paper label sets (tool-independent; needs only downloaded data).
#SBATCH --clusters=cbsubscb
#SBATCH --account=danko_0001
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=ppb_labels
#SBATCH --output=/home/ayh8/slurm/ppb_labels.%j.out

source /home2/ayh8/pydreg_pints/slurm/common.sh
make_scratch
cd "$PROJ"
"$VENV/python" scripts/01_build_labels.py
