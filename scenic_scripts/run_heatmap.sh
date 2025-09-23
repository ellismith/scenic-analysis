#!/usr/bin/env bash
set -euo pipefail

PYBIN="/scratch/easmit31/conda_envs/pyscenic/bin/python"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <auc_mtx.csv>"
    exit 1
fi

INPUT="$1"
BASENAME=$(basename "$INPUT" .csv)   # e.g. astros_adults_allLv_auc_mtx
PREFIX=${BASENAME%_auc_mtx}          # strip suffix -> astros_adults_allLv

echo "[heatmap] Input:   $INPUT"
echo "[heatmap] Prefix:  $PREFIX"

# Count dimensions
CELLS=$(tail -n +2 "$INPUT" | wc -l)
REGS=$(awk -F, 'NR==1{print NF-1}' "$INPUT")
echo "[heatmap] AUCell matrix: ${CELLS} cells × ${REGS} regulons"

# Summary
$PYBIN plot_aucell_summary.py \
  --input "$INPUT" \
  --out "${PREFIX}_regulon_summary.csv"

# Heatmap
$PYBIN plot_aucell_heatmap.py \
  --input "$INPUT" \
  --out "${PREFIX}_auc_heatmap.png" \
  --n-cells 200 \
  --n-regs 30 \
  --seed 0

echo "[heatmap] Artifacts written:"
ls -lh "${PREFIX}_regulon_summary.csv" "${PREFIX}_auc_heatmap.png"
