#!/bin/bash
#SBATCH --job-name=expr_umap
#SBATCH --output=logs/expr_umap_%j.out
#SBATCH --error=logs/expr_umap_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

# Load conda
source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
conda activate /scratch/easmit31/conda_envs/pyscenic

PYBIN=/scratch/easmit31/conda_envs/pyscenic/bin/python
PROJECT_DIR=/scratch/easmit31/GRN_copy/scenic/scenic_scripts

cd $PROJECT_DIR

$PYBIN 85_visualize_expr_umap.py \
  --input /scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_adults_allLv.loom \
  --colorby age \
  --out /scratch/easmit31/GRN_copy/scenic/figures/astros_expr_umap.png
