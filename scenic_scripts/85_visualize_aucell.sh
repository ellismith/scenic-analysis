#!/bin/bash
#SBATCH --job-name=viz
#SBATCH --output=/home/easmit31/scenic/logs/85_visualize_aucell_%j.out
#SBATCH --error=/home/easmit31/scenic/logs/85_visualize_aucell_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
set -euo pipefail
set -x

PY="/scratch/easmit31/conda_envs/pyscenic/bin/python"
cd /home/easmit31/scenic

# Forward all command-line args to Python script
$PY 85_visualize_aucell.py "$@"
