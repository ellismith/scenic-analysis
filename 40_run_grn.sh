#!/bin/bash
#SBATCH --job-name=grn
#SBATCH --time=03:59:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/40_run_grn_%j.out
#SBATCH --partition=htc

# load mamba and activate your conda environment
module load gcc-11.2.0-gcc-8.5.0
module load mamba/latest
source activate /scratch/easmit31/conda_envs/pyscenic_final
PROJECT_DIR="/scratch/easmit31/GRN/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic_final/bin/python"
cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"
mkdir -p logs
# run the step, forwarding all args
"$PYTHON" 40_run_grn.py "$@"
