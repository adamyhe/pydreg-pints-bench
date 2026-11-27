#!/bin/bash -l
# Verbose report from bootstrap CIs (run after 40_bootstrap).
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=ppb_report

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
cd "$PROJ"
"$VENV/python" scripts/04_summarize.py > /dev/null
"$VENV/python" scripts/06_report.py
