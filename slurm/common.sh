# Sourced by the sbatch scripts. Resolves the project, venv and the array task's library.
set -euo pipefail
PROJ=/home2/ayh8/pydreg_pints
VENV="$PROJ/.venv/bin"
export PPBENCH_PROJ="$PROJ"
export PYTHONNOUSERSITE=1

# lib for this array index (1-based over the rows of config/libraries.tsv in phase $PPB_PHASE:
# 1 (default), 2, or all)
PPB_PHASE="${PPB_PHASE:-1}"
lib_for_task() {
    awk -F'\t' -v i="$1" -v ph="$PPB_PHASE" \
        'NR>1 && (ph=="all" || $2==ph) {n++; if (n==i) {print $1; exit}}' "$PROJ/config/libraries.tsv"
}
lib_field() {  # lib_field LIB COLUMN_NAME
    awk -F'\t' -v lib="$1" -v col="$2" 'NR==1 {for (i=1;i<=NF;i++) if ($i==col) c=i; next} $1==lib {print $c}' \
        "$PROJ/config/libraries.tsv"
}

make_scratch() {
    SCRATCH=/workdir/$USER/$SLURM_JOB_ID
    mkdir -p "$SCRATCH"
    export TMPDIR="$SCRATCH/tmp"
    mkdir -p "$TMPDIR"
    trap 'rm -rf "$SCRATCH"' EXIT
}
