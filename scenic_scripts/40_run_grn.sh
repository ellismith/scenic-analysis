#!/bin/bash
#SBATCH --job-name=grn
#SBATCH --time=1-00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/40_run_grn_%j.out
#SBATCH --partition=highmem

# load mamba and activate your conda environment
module load mamba/latest
source activate /scratch/easmit31/conda_envs/pyscenic_final

PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic_final/bin/python"

cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"
mkdir -p logs

# run the step, forwarding all args
"$PYTHON" 40_run_grn.py "$@"
