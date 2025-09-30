#!/bin/bash
#SBATCH --job-name=sel
#SBATCH --time=0-10
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/easmit31/GRN_copy/scenic/scenic_scripts/logs/10_select_cells_%j.out
#SBATCH --error=/scratch/easmit31/GRN_copy/scenic/scenic_scripts/logs/10_select_cells_%j.err

PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
LOG_DIR="$PROJECT_DIR/logs"

/bin/mkdir -p "$LOG_DIR"

cd "$PROJECT_DIR" || { echo "Project dir not found: $PROJECT_DIR"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

echo "[\$(date)] Running: /scratch/easmit31/conda_envs/pyscenic_final/bin/python 10_select_cells.py \$@"
/scratch/easmit31/conda_envs/pyscenic_final/bin/python 10_select_cells.py "$@"
