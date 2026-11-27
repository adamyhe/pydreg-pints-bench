#!/bin/bash -l
# Memory profile of pydreg (as benchmarked: cupy backend, 16 processes) on K562_PROseq_G1 and K562_PROseq_G7,
# for comparison with 50_memprof_pints. Calls are discarded; only results/memprof/ is kept.
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tP100:1
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=1-2%2
#SBATCH --job-name=ppb_memprof_pydreg

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
LIBS=(K562_PROseq_G1 K562_PROseq_G7); LIB=${LIBS[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "lib=$LIB host=$(hostname) gpu=$CUDA_VISIBLE_DEVICES"

cp "$PROJ/data/libs/$LIB.pl.bw" "$PROJ/data/libs/$LIB.mn.bw" "$SCRATCH/"
mkdir -p "$SCRATCH/out" "$PROJ/results/memprof"
cd "$SCRATCH"
"$VENV/python" "$PROJ/scripts/08_memprof.py" "$PROJ/results/memprof/$LIB.pydreg.p16" -- \
    "$VENV/pydreg" --backend cupy -p "$SLURM_NTASKS" --no-progress "$LIB.pl.bw" "$LIB.mn.bw" "out/$LIB"
