#!/bin/bash
#SBATCH --job-name=age_sex_scatter
#SBATCH --output=logs/73_scatter_%j.out
#SBATCH --error=logs/73_scatter_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
set -euo pipefail
set -x

# Load mamba
module load mamba/latest

# Define project + python
PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic/bin/python"

cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Args: LM_TSV_PATH FIGURES_DIR
"$PYTHON" 73_age_sex_scatter.py "$1" "$2"
