# Site settings, sourced by the fetch scripts, slurm/common.sh and slurm/sub.sh.
# Loads config/site.env if present (copy config/site.env.example), then fills defaults.
PPB_PROJ="${PPB_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
if [[ -f "$PPB_PROJ/config/site.env" ]]; then
    set -a; source "$PPB_PROJ/config/site.env"; set +a
fi
export PPB_PROJ PPBENCH_PROJ="$PPB_PROJ"
export PPB_SLURM_LOGS="${PPB_SLURM_LOGS:-$PPB_PROJ/logs/slurm}"
export PPB_SCRATCH_ROOT="${PPB_SCRATCH_ROOT:-/workdir/$USER}"
if [[ -z "${KENT_BIN:-}" ]]; then
    if command -v bedGraphToBigWig >/dev/null 2>&1; then KENT_BIN="$(dirname "$(command -v bedGraphToBigWig)")"
    else KENT_BIN=/programs/kentUtils/bin; fi
fi
export KENT_BIN
export PPB_SLURM_CLUSTER="${PPB_SLURM_CLUSTER:-}" PPB_SLURM_ACCOUNT="${PPB_SLURM_ACCOUNT:-}"
export PPB_SBATCH_ARGS="${PPB_SBATCH_ARGS:-}" PPB_GPU_SBATCH_ARGS="${PPB_GPU_SBATCH_ARGS:-}"
export DREG_K562_HMM_BED="${DREG_K562_HMM_BED:-}" GM12878_UW_DHS_BED="${GM12878_UW_DHS_BED:-}"
