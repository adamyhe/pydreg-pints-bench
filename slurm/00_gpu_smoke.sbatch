#!/bin/bash -l
# Quick check that CuPy sees the P100 and pydreg selects the cupy backend. Run this first.
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tP100:1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --job-name=ppb_gpu_smoke

source "${PPB_PROJ:-$SLURM_SUBMIT_DIR}/slurm/common.sh"   # submit with slurm/sub.sh
nvidia-smi
"$VENV/python" - <<'EOF'
import cupy as cp
print("cupy", cp.__version__, "devices", cp.cuda.runtime.getDeviceCount(),
      "cc", cp.cuda.Device(0).compute_capability)
a = cp.random.rand(2048, 2048, dtype=cp.float32)
print("matmul ok", float((a @ a).sum()) > 0)
from pydreg import backend
print("pydreg auto-detected backend:", backend.detect_backend())
EOF
"$VENV/pydreg" --help | head -5
