#!/bin/bash
#SBATCH --job-name=aucell_barplots
#SBATCH --output=logs/85_barplots_%j.out
#SBATCH --error=logs/85_barplots_%j.err
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

# Run barplot script with args
"$PYTHON" 85_visualize_aucell.py \
  --input "$1" \
  --figdir "$2"
