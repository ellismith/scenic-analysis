#!/bin/bash
#SBATCH --job-name=subset_adults
#SBATCH --output=logs/25_subset_adults_%j.out
#SBATCH --error=logs/25_subset_adults_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=highmem

set -euo pipefail
set -x

# Load conda without sourcing bashrc
source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
conda activate /scratch/easmit31/conda_envs/pyscenic

PYTHON=/scratch/easmit31/conda_envs/pyscenic/bin/python

$PYTHON 25_subset_adults.py "$@"
