#!/bin/bash -l
# Verbose report from bootstrap CIs (run after 40_bootstrap).
#SBATCH --clusters=cbsubscb
#SBATCH --account=danko_0001
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=ppb_report
#SBATCH --output=/home/ayh8/slurm/ppb_report.%j.out

source /home2/ayh8/pydreg_pints/slurm/common.sh
cd "$PROJ"
"$VENV/python" scripts/04_summarize.py > /dev/null
"$VENV/python" scripts/06_report.py
