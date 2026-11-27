#!/bin/bash -l
# pydreg on each phase-1 library (one GPU per task; at most 2 tasks at once = the node's 2 P100s).
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tP100:1
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-5%2
#SBATCH --job-name=ppb_pydreg

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
LIB=$(lib_for_task "$SLURM_ARRAY_TASK_ID")
echo "lib=$LIB host=$(hostname) gpu=$CUDA_VISIBLE_DEVICES"
nvidia-smi --query-gpu=name,memory.total --format=csv

cp "$PROJ/data/libs/$LIB.pl.bw" "$PROJ/data/libs/$LIB.mn.bw" "$SCRATCH/"
mkdir -p "$SCRATCH/out"
cd "$SCRATCH"
/usr/bin/time -v "$VENV/pydreg" --backend cupy -p "$SLURM_NTASKS" -v --no-progress \
    "$LIB.pl.bw" "$LIB.mn.bw" "out/$LIB" 2> >(tee "out/$LIB.pydreg.stderr.log" >&2)

DEST="$PROJ/results/calls/$LIB/pydreg"
mkdir -p "$DEST"
cp out/* "$DEST/"
echo "copied to $DEST"; ls -la "$DEST"
