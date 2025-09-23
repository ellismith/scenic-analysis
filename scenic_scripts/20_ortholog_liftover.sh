#!/bin/bash
#SBATCH --job-name=liftover
#SBATCH --partition=highmem
#SBATCH --time=24:00:00
#SBATCH --mem=256GB
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/20_ortholog_liftover_%j.out

# load mamba
module load mamba/latest

PROJECT_DIR="/scratch/easmit31/GRN_copy/scenic/scenic_scripts"
PYTHON="/scratch/easmit31/conda_envs/pyscenic/bin/python"

cd "$PROJECT_DIR" || { echo "Project dir not found"; exit 1; }
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"

# Choose which Python script to run:
SCRIPT="20_ortholog_liftover.py"   # default

if [[ "$1" == "--script" ]]; then
    shift
    if [[ "$1" == "sparse" ]]; then
        SCRIPT="20_ortholog_liftover_sparse.py"
    elif [[ "$1" == "dense" ]]; then
        SCRIPT="20_ortholog_liftover.py"
    fi
    shift
fi

echo "Running liftover with: $SCRIPT"

if [[ "$SCRIPT" == "20_ortholog_liftover_sparse.py" ]]; then
    # Pass only the 3 arguments (skip the script name duplication)
    "$PYTHON" "$SCRIPT" "$@"
else
    "$PYTHON" "$SCRIPT" "$@"
fi
