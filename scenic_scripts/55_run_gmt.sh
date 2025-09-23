#!/bin/bash
#SBATCH --job-name=gmt
#SBATCH --output=logs/55_run_gmt_%j.out
#SBATCH --error=logs/55_run_gmt_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail
set -x

# pySCENIC executable
export PATH="/scratch/easmit31/conda_envs/pyscenic/bin:$PATH"
PYS=/scratch/easmit31/conda_envs/pyscenic/bin/pyscenic
PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Convert regulons CSV to GMT format (collapses duplicate TF regulons)
$PYS regulons-to-gmt \
  /scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_allLv_reg.csv \
  --output /scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_allLv_reg.gmt
