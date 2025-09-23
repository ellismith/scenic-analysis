#!/bin/bash
#SBATCH --job-name=plots95
#SBATCH --time=0-04:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/95_plots_%j.out

module load mamba/latest
# activate your conda environment
source activate /scratch/easmit31/conda_envs/pyscenic

# absolute paths
PROJECT_DIR="/home/easmit31/scenic"
PYTHON="/scratch/easmit31/conda_envs/pyscenic/bin/python"
LOG_DIR="/home/easmit31/scenic/logs"

# 1) go to project and ensure logs dir exists there
cd "$PROJECT_DIR" || { echo "Project dir not found: $PROJECT_DIR"; exit 1; }
mkdir -p "$LOG_DIR"

# 2) make sure Python can import config.py from this folder
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

'$PYTHON' 95_plots.py --top_n 5
