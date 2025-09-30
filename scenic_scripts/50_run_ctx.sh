#!/bin/bash
#SBATCH --job-name=ctx
#SBATCH --output=logs/50_run_ctx_%j.out
#SBATCH --error=logs/50_run_ctx_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

set -euo pipefail
set -x

# pySCENIC environment
PYTHON=/scratch/easmit31/conda_envs/pyscenic_final/bin/python
PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"

cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Call the Python wrapper which handles glob expansion
"$PYTHON" 50_run_ctx.py "$@"
