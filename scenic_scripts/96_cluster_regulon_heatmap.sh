#!/bin/bash
#SBATCH --job-name=cluster_heatmap
#SBATCH --output=logs/96_cluster_heatmap_%j.out
#SBATCH --error=logs/96_cluster_heatmap_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
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

# Args: INPUT CLUSTER_ID CELLTYPE FIGDIR MAX_REGULONS
"$PYTHON" 96_cluster_regulon_heatmap.py \
  --input "$1" \
  --cluster "$2" \
  --celltype "$3" \
  --figdir "$4" \
  --max_regulons "$5"
