#!/bin/bash -l
# Memory profile of PINTS (scripts/08_memprof.py: PSS summed over the process tree, sampled every second).
# Tasks: {K562_PROseq_G1 (deepest), K562_PROseq_G7} x {--debug 16 threads (as benchmarked), no --debug 16 threads,
# --debug 4 threads}. Calls are discarded; only results/memprof/ is kept.
#SBATCH --partition=regular,long7,long30
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --array=1-6
#SBATCH --job-name=ppb_memprof_pints

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
make_scratch
LIBS=(K562_PROseq_G1 K562_PROseq_G7); CFGS=(debug_t16 nodebug_t16 debug_t4)
i=$((SLURM_ARRAY_TASK_ID - 1)); LIB=${LIBS[$((i / 3))]}; CFG=${CFGS[$((i % 3))]}
EXP=$(lib_field "$LIB" pints_exp_type)
THREADS=${CFG##*_t}; DEBUG=$([[ $CFG == debug_* ]] && echo --debug || true)
echo "lib=$LIB cfg=$CFG exp_type=$EXP host=$(hostname)"

cp "$PROJ/data/libs/$LIB.pl.bw" "$PROJ/data/libs/$LIB.mn.bw" "$SCRATCH/"
mkdir -p "$SCRATCH/out" "$PROJ/results/memprof"
cd "$SCRATCH"
"$VENV/python" "$PROJ/scripts/08_memprof.py" "$PROJ/results/memprof/$LIB.pints.$CFG" -- \
    "$VENV/pints_caller" --bw-pl "$LIB.pl.bw" --bw-mn "$LIB.mn.bw" --exp-type "$EXP" \
    --save-to out --file-prefix "$LIB" --thread "$THREADS" \
    --chromosome-start-with chr -f chrM _ EBV $DEBUG --dont-check-updates
