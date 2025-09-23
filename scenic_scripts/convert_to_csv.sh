#!/usr/bin/env bash
set -euo pipefail

# Convert AUCell .tsv scores into .csv so plotting scripts can read them.

for TSV in /scratch/easmit31/GRN_copy/scenic/h5ad_files/*_aucell_scores.tsv; do
    if [[ -f "$TSV" ]]; then
        PREFIX=$(basename "$TSV" _aucell_scores.tsv)
        CSV="${PREFIX}_auc_mtx.csv"
        echo "[convert] Converting $TSV -> $CSV"
        # Replace tabs with commas
        sed 's/\t/,/g' "$TSV" > "/scratch/easmit31/GRN_copy/scenic/h5ad_files/$CSV"
    fi
done

echo "[convert] Done. Available CSVs:"
ls -lh /scratch/easmit31/GRN_copy/scenic/h5ad_files/*_auc_mtx.csv
