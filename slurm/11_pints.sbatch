#!/bin/bash -l
# PINTS on each phase-1 library. --debug writes *_peaks.gtf with every per-strand candidate and its
# p-value (needed to sweep PINTS over its full range, as the PINTS paper's own ROC did).
#SBATCH --partition=regular,long7,long30
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=1-5
#SBATCH --job-name=ppb_pints

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
LIB=$(lib_for_task "$SLURM_ARRAY_TASK_ID")
EXP=$(lib_field "$LIB" pints_exp_type)
echo "lib=$LIB exp_type=$EXP host=$(hostname)"

cp "$PROJ/data/libs/$LIB.pl.bw" "$PROJ/data/libs/$LIB.mn.bw" "$SCRATCH/"
mkdir -p "$SCRATCH/out"
cd "$SCRATCH"
/usr/bin/time -v "$VENV/pints_caller" \
    --bw-pl "$LIB.pl.bw" --bw-mn "$LIB.mn.bw" --exp-type "$EXP" \
    --save-to out --file-prefix "$LIB" --thread "$SLURM_NTASKS" \
    --chromosome-start-with chr -f chrM _ EBV \
    --debug --dont-check-updates 2> >(tee "out/$LIB.pints.stderr.log" >&2)

DEST="$PROJ/results/calls/$LIB/pints"
mkdir -p "$DEST"
cp -r out/* "$DEST/"
echo "copied to $DEST"; ls -la "$DEST"
