#!/bin/bash
#SBATCH --job-name=umap
#SBATCH --output=logs/plot_aucell_umap_%j.out
#SBATCH --error=logs/plot_aucell_umap_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
set -euo pipefail
set -x

PY="/scratch/easmit31/conda_envs/pyscenic/bin/python"
cd /scratch/easmit31/GRN_copy/scenic/scenic_scripts

# Forward all command-line args directly to the Python script
$PY plot_aucell_umap.py "$@"
