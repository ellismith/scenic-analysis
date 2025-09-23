#!/bin/bash
#SBATCH --job-name=allclusters_heatmaps
#SBATCH --output=logs/96_allclusters_%j.out
#SBATCH --error=logs/96_allclusters_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
set -euo pipefail
set -x

# Load mamba
module load mamba/latest

# Define project + python
PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic/bin/python"

cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Run the script with args
"$PYTHON" 96_all_clusters_regulon_heatmaps.py \
  --input "$1" \
  --celltype "$2" \
  --figdir "$3" \
  --max_regulons "$4"
