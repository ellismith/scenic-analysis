#!/usr/bin/env python3
"""
00_config.py
Central configuration for the GRN/SCENIC pipeline.
Keeps all file paths, dataset identifiers, and parameters in one place
so they are consistent across scripts.
Usage:
    from config import H5AD_IN, WORK_DIR, N_JOBS, ...
"""

from pathlib import Path
import os

# -------------------
# Core labels
# -------------------
# Cell type / dataset tag (can override with env var CT)
CT = os.environ.get("CT", "glutamatergic-neurons_L2-3_TEST")

# Louvain clusters of interest (override with env var SELECTED_LV="1,2,3")
_selected_lv = os.environ.get("SELECTED_LV", "")
if _selected_lv:
    SELECTED_LV = [int(x) for x in _selected_lv.split(",")]
else:
    SELECTED_LV = [5, 0, 17, 27, 29, 16, 24, 30, 28]

# -------------------
# Base directories
# -------------------
PROCESSED_DIR = Path(
    os.environ.get("PROCESSED_DIR", "/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update")
)

WORK_DIR = Path(os.environ.get("WORK_DIR", f"/scratch/easmit31/GRN_copy/{CT}"))
WORK_DIR.mkdir(parents=True, exist_ok=True)

# -------------------
# Input / mapping files
# -------------------
# H5AD input: must be set explicitly since filenames vary
H5AD_IN = Path(
    os.environ.get(
        "H5AD_IN",
        str(PROCESSED_DIR / f"{CT}.h5ad"),  # default guess, but override recommended
    )
)

# Ortholog mapping file
MART_CSV = Path(
    os.environ.get(
        "MART_CSV",
        "/scratch/wyang114/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv",
    )
)

# -------------------
# Outputs
# -------------------
SUBSET_H5AD   = WORK_DIR / f"{CT}_subset.h5ad"
LIFTOVER_H5AD = WORK_DIR / f"{CT}_liftover_adata.h5ad"
SCENIC_LOOM   = WORK_DIR / f"{CT}_filtered_scenic.loom"
ADJ_CSV       = WORK_DIR / "adj.csv"
REG_CSV       = WORK_DIR / f"{CT}_reg.csv"
PYSCENIC_OUT  = WORK_DIR / f"{CT}_pyscenic_output.h5ad"

# -------------------
# PySCENIC resources
# -------------------
TF_LIST     = Path("/scratch/wyang114/data/cisTarget_databases/tf_lists/allTFs_hg38.txt")
RANK_DB_GLOB = os.environ.get(
    "RANK_DB_GLOB",
    "/scratch/wyang114/data/cisTarget_databases/mc_v10_clust_ranking/human/*.feather",
)
MOTIF_ANNOT = Path(
    os.environ.get(
        "MOTIF_ANNOT",
        "/scratch/wyang114/data/cisTarget_databases/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    )
)

# -------------------
# Performance / filtering parameters
# -------------------
NCELLS_FILTER = int(os.environ.get("NCELLS_FILTER", 100))
N_JOBS        = int(os.environ.get("N_JOBS", 20))

# Logs folder for Slurm outputs
LOG_DIR = WORK_DIR / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
