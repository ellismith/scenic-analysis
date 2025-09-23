#!/bin/bash
#SBATCH --job-name=ctx
#SBATCH --output=logs/50_run_ctx_%j.out
#SBATCH --error=logs/50_run_ctx_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

set -euo pipefail
set -x

# pySCENIC environment
export PATH="/scratch/easmit31/conda_envs/pyscenic/bin:$PATH"
PYS=/scratch/easmit31/conda_envs/pyscenic/bin/pyscenic
PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Run ctx, forwarding all arguments
$PYS ctx "$@"
