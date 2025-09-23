#!/bin/bash
#SBATCH --job-name=viz_age
#SBATCH --output=logs/viz_age_%j.out
#SBATCH --error=logs/viz_age_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

# Load conda safely (avoid sourcing ~/.bashrc to skip /etc/bashrc issues)
source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
conda activate /scratch/easmit31/conda_envs/pyscenic

PYBIN=/scratch/easmit31/conda_envs/pyscenic/bin/python
PROJECT_DIR=/scratch/easmit31/GRN_copy/scenic/scenic_scripts

cd $PROJECT_DIR

# Forward args directly to the Python script
$PYBIN "$@"
