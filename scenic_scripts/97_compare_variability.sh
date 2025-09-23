#!/bin/bash
#SBATCH --job-name=compare_variability
#SBATCH --output=logs/97_compare_variability_%j.out
#SBATCH --error=logs/97_compare_variability_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

# Load mamba
module load mamba/latest

# Define project + python
PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic/bin/python"

cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Args: H5AD1 H5AD2 --topN N --out OUT
"$PYTHON" 97_compare_variability.py "$@"
