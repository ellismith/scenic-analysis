#!/bin/bash
#SBATCH --job-name=aucell_full
#SBATCH --output=aucell_full_%j.out
#SBATCH --error=aucell_full_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G

set -euo pipefail
export PYTHONNOUSERSITE=1

# ===================== USER PARAMS =====================
# Path to the big h5ad
H5AD="/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/Res1_glutamatergic-neurons_update.h5ad"

# Regulon definitions you already used (CSV with TF,target OR a .gmt)
REGULONS="regulons_for_aucell.csv"   # or: regulons_for_aucell.gmt

# Subsample size: set to 0 to use ALL CELLS; set to e.g. 100000 to downsample
N_CELLS=0

# Output prefix (files land in current folder)
PREFIX="aucell_full"
# ======================================================

PYBIN="/scratch/easmit31/conda_envs/pyscenic/bin/python"
PYS="/scratch/easmit31/conda_envs/pyscenic/bin/pyscenic"

echo "[INFO] Job started on $(hostname) at $(date)"
echo "[INFO] H5AD:      $H5AD"
echo "[INFO] REGULONS:  $REGULONS"
echo "[INFO] N_CELLS:   $N_CELLS (0 means ALL)"
echo "[INFO] PREFIX:    $PREFIX"

# 1) Make a loom from the .h5ad (all cells or subsample)
#    We do this in Python so it’s reproducible and logs the actual counts used.
LOOM="${PREFIX}.loom"

$PYBIN - <<PY
import scanpy as sc, numpy as np, os, sys

h5ad = "${H5AD}"
out_loom = "${LOOM}"
n_cells = int("${N_CELLS}")

print("[PY] Reading AnnData (backed='r'):", h5ad, flush=True)
adata = sc.read_h5ad(h5ad, backed="r")
n_obs, n_vars = adata.n_obs, adata.n_vars
print(f"[PY] AnnData dims: cells={n_obs}, genes={n_vars}", flush=True)

# Decide which cells to keep
if n_cells and n_cells < n_obs:
    # Random subsample
    rng = np.random.default_rng(0)
    keep_idx = rng.choice(n_obs, size=n_cells, replace=False)
    # Need to materialize a view for write (this will load the slice in RAM)
    # With 300G this should be okay for ~100k cells; for ALL, we skip subsample.
    print(f"[PY] Subsampling to {n_cells} cells.", flush=True)
    adata = sc.read_h5ad(h5ad)  # load into memory once for slicing/write
    adata = adata[keep_idx, :].copy()
else:
    print("[PY] Using ALL cells.", flush=True)
    # For write_loom we need an in-memory object; backed views generally can't be written.
    adata = sc.read_h5ad(h5ad)  # load into memory (we requested 300G)

print("[PY] Writing loom:", out_loom, flush=True)
adata.write_loom(out_loom)
print("[PY] Done loom write.", flush=True)
PY

# 2) Run AUCell on the loom
#    Output is a CSV of AUCell scores: rows=cells, cols=regulons
AUC_CSV="${PREFIX}_auc_mtx.csv"

echo "[INFO] Running AUCell…"
# Note: remove --verbose (not supported in your pyscenic version)
$PYS aucell "${LOOM}" "${REGULONS}" -o "${AUC_CSV}" --num_workers 32

echo "[INFO] Done AUCell. Files written:"
ls -lh "${LOOM}" "${AUC_CSV}" 2>/dev/null || true
echo "[INFO] Finished at $(date)"
