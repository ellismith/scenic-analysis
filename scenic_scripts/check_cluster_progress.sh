#!/bin/bash
# Usage: ./check_cluster_progress.sh <celltype> <nclusters>
# Example: ./check_cluster_progress.sh gaba 26

if [ $# -lt 2 ]; then
  echo "❌ Usage: $0 <celltype> <nclusters>"
  exit 1
fi

CELLTYPE=$1
NCLUSTERS=$2
FIG_DIR="/scratch/easmit31/GRN_copy/scenic/figures"
LOG_DIR="logs"

for cl in $(seq 0 $((NCLUSTERS-1))); do
  PNG="${FIG_DIR}/${CELLTYPE}_cluster${cl}_regulons.png"
  OUT_LOG=$(ls -1 ${LOG_DIR}/${CELLTYPE}_cluster${cl}_heatmap_*.out 2>/dev/null | tail -n1)
  ERR_LOG=$(ls -1 ${LOG_DIR}/${CELLTYPE}_cluster${cl}_heatmap_*.err 2>/dev/null | tail -n1)

  if [ -f "$PNG" ]; then
    echo "cluster ${cl} → ✅ PNG exists: $PNG"
  elif squeue -u $USER | grep -q "${CELLTYPE}_cl${cl}"; then
    echo "cluster ${cl} → ⏳ still running"
  elif [ -f "$ERR_LOG" ] && grep -qi "error" "$ERR_LOG"; then
    echo "cluster ${cl} → ❌ failed (see $ERR_LOG)"
  else
    echo "cluster ${cl} → ❓ no PNG, not running (check logs)"
  fi
done
