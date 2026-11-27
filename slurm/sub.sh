#!/bin/bash
# Submit one sbatch script with the site settings; prints the job ID.
#   bash slurm/sub.sh [sbatch options] slurm/<script>.sbatch
# Adds --clusters/--account (if set), --output under $PPB_SLURM_LOGS, then $PPB_SBATCH_ARGS and, for GPU
# scripts, $PPB_GPU_SBATCH_ARGS; your own options come last and win.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/scripts/site.sh"
script="${!#}"
mkdir -p "$PPB_SLURM_LOGS"
args=()
[[ -n "$PPB_SLURM_CLUSTER" ]] && args+=(--clusters="$PPB_SLURM_CLUSTER")
[[ -n "$PPB_SLURM_ACCOUNT" ]] && args+=(--account="$PPB_SLURM_ACCOUNT")
if grep -q -- "--array" "$script" || [[ " $* " == *" --array"* ]]; then
    args+=(--output="$PPB_SLURM_LOGS/%x.%A_%a.out")
else
    args+=(--output="$PPB_SLURM_LOGS/%x.%j.out")
fi
read -r -a extra <<< "$PPB_SBATCH_ARGS"; args+=("${extra[@]}")
if grep -q -- "--gres=gpu" "$script"; then read -r -a extra <<< "$PPB_GPU_SBATCH_ARGS"; args+=("${extra[@]}"); fi
cd "$PPB_PROJ"
sbatch --parsable --export=ALL "${args[@]}" "$@" | cut -d';' -f1
