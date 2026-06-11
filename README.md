# SCENIC Gene Regulatory Network Analysis Pipeline
## Snyder Lab | ASU | Macaque Brain Aging
 
Pipeline for inferring transcription factor gene regulatory networks (GRNs) from single-nucleus RNA-seq data using [pySCENIC](https://pyscenic.readthedocs.io/), applied to adult rhesus macaque brain across 11 regions and multiple cell types.
 
---
 
## Overview
 
This pipeline:
- Lifts over macaque gene expression to human orthologs (required for cisTarget databases)
- Infers TF co-expression networks using GRNBoost2
- Prunes networks using cisTarget motif databases
- Scores per-cell regulon activity using AUCell
- Models regulon activity as a function of age and sex using mixed-effects models, with animal as a random intercept
---
 
## Repository Structure
 
```
scenic-analysis/
├── scenic_scripts/              # All pipeline scripts
│   ├── 10_select_cells.*        # Filter to adults + cell type
│   ├── 20_ortholog_liftover.*   # Macaque → human gene space
│   ├── 20_ortholog_liftover_sparse.py  # Canonical liftover script
│   ├── 30_convert_to_loom.*     # h5ad → loom for GRNBoost
│   ├── 40_run_grn.*             # GRNBoost2 network inference
│   ├── 50_run_ctx.*             # cisTarget motif pruning
│   ├── 56_convert_regulons.py   # Parse regulons → target pairs
│   ├── 60_aucell_final.*        # AUCell scoring
│   ├── 70_regulon_mixed_models.*   # Mixed-effects models (global)
│   ├── 74_cluster_mixed_models.*   # Mixed-effects models (per cluster)
│   ├── chunk_h5ad.py            # Chunk large h5ads for AUCell
│   ├── merge_aucell_chunks.*    # Merge chunked AUCell outputs
│   ├── plotting_scripts/        # Visualization scripts
│   ├── not_using/               # Archived old scripts
│   └── references/              # TF lists, regulon summaries
```
 
On Sol, all data lives under `/scratch/easmit31/GRN/scenic/`:
```
scenic/
├── h5ad_adult_files/        # Per-cell-type h5ads
├── adult_grn_outputs/       # GRN adjacencies and regulon CSVs
└── results/                 # Mixed model outputs
```
 
---
 
## Resources (Sol HPC paths)
 
| Resource | Path |
|----------|------|
| TF list (hg38) | `/scratch/easmit31/GRN/scenic/allTFs_hg38.txt` |
| Ranking DB (10kb) | `/scratch/easmit31/data/cisTarget_databases/mc_v10_clust_ranking/human/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` |
| Ranking DB (500bp) | `/scratch/easmit31/data/cisTarget_databases/mc_v10_clust_ranking/human/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` |
| Motif annotations | `/scratch/easmit31/data/cisTarget_databases/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` |
| Conda environment | `/scratch/easmit31/conda_envs/pyscenic_final` |
 
---
 
## Pipeline Steps
 
### Step 10 — Select cells
Filter to a single cell type, adults only (age ≥ 1 year).
```bash
sbatch 10_select_cells.sh [args]
```
 
### Step 20 — Ortholog liftover
Map macaque Ensembl IDs to human HGNC symbols (1:1 orthologs only). Optionally downsample cells.
```bash
sbatch 20_ortholog_liftover.sh \
  <input_selected.h5ad> \
  <mapping.csv> \
  <output_liftover.h5ad>
```
To use all cells (no downsampling): verify `--downsample` behavior in `20_ortholog_liftover_sparse.py`.
 
### Step 30 — Convert to loom
Convert h5ad to loom format required by GRNBoost2.
```bash
sbatch 30_convert_to_loom.sh \
  --adata_in <liftover.h5ad> \
  --loom_out <output.loom>
```
 
### Step 40 — GRN inference (GRNBoost2)
Infer TF-gene co-expression network.
```bash
sbatch 40_run_grn.sh \
  --loom_in <input.loom> \
  --tfs /scratch/easmit31/GRN/scenic/allTFs_hg38.txt \
  --out <output_adj.csv>
```
 
### Step 50 — cisTarget pruning
Prune adjacencies using cisTarget motif databases to get regulons.
```bash
sbatch 50_run_ctx.sh \
  --adj <adj.csv> \
  --expression_mtx_fname <input.loom> \
  --db_glob "/scratch/easmit31/data/cisTarget_databases/mc_v10_clust_ranking/human/hg38_*.rankings.feather" \
  --annotations_fname /scratch/easmit31/data/cisTarget_databases/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --output <regulons.csv>
```
 
### Step 56 — Convert regulons
Parse raw cisTarget CSV into simple `regulon,target` pairs for AUCell.
```bash
/scratch/easmit31/conda_envs/pyscenic_final/bin/python 56_convert_regulons.py \
  <regulons.csv> \
  <regulons_simple.csv>
```
 
### Step 60 — AUCell scoring
Score each cell for regulon activity. Uses Python API directly (not CLI) to avoid sparse matrix issues.
```bash
sbatch 60_aucell_final.sh \
  <liftover.h5ad> \
  <regulons_simple.csv> \
  <aucell_output.h5ad>
```
For large cell types (>200k cells), chunk first using `chunk_h5ad.py` + `chunk_job.sh`, then merge with `merge_aucell_chunks.sh`.
 
### Step 70 — Global mixed-effects models
Model regulon activity as a function of age and sex across all cells.
```bash
sbatch 70_regulon_mixed_models.sh \
  <aucell.h5ad> \
  <output_mixed_models.tsv>
```
Formula: `activity_std ~ age_std + sex` with random intercept per `animal_id`. Age and activity are z-scored. FDR correction (BH) applied per term.
 
### Step 74 — Cluster-level mixed-effects models
Same models run per louvain cluster. Sex-chromosome TFs excluded using chromosome annotations from source h5ad.
```bash
sbatch 74_cluster_mixed_models.sh \
  <aucell.h5ad> \
  <selected.h5ad> \
  <output_cluster_mixed_models.tsv>
```
 
---
 
## File Naming Convention
 
```
adult_grn_<celltype>_selected.h5ad           # Step 10 output
adult_grn_<celltype>_liftover.h5ad           # Step 20 output
adult_grn_<celltype>.loom                    # Step 30 output
adult_grn_<celltype>_adj.csv                 # Step 40 output
adult_grn_<celltype>_regulons.csv            # Step 50 output
adult_grn_<celltype>_regulons_simple.csv     # Step 56 output
adult_grn_<celltype>_aucell.h5ad             # Step 60 output
adult_grn_<celltype>_mixed_models.tsv        # Step 70 output
adult_grn_<celltype>_cluster_mixed_models.tsv  # Step 74 output
```
 
---
 
## Cell Type Status
 
| Cell type | Cells | Select | Liftover | GRN | CTX | AUCell | Mixed models |
|-----------|-------|--------|----------|-----|-----|--------|--------------|
| Astrocytes | ~416k | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| GABAergic | ~523k | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Microglia | ~158k (100k liftover) | ✓ | ✓ | 🔄 | — | — | — |
| MSN | ~157k | ✓ | ✓ | — | — | — | — |
| OPC/Oligo | ~1.04M | ✓ | ✓ | — | — | — | — |
| Glutamatergic | TBD | ✓ | — | — | — | — | — |
| Basket cells | — | skipped | | | | | |
 
---
 
## Environment Setup
 
**Environment**: `/scratch/easmit31/conda_envs/pyscenic_final`
 
**Critical**: load gcc module before any pySCENIC work:
```bash
module load gcc-11.2.0-gcc-8.5.0
```
 
**Rebuild from scratch**:
```bash
module load mamba/latest
mamba env create -p /scratch/easmit31/conda_envs/pyscenic_final \
  -f /scratch/easmit31/conda_envs/pyscenic_env.yaml
# Then fix setuptools:
/scratch/easmit31/conda_envs/pyscenic_final/bin/pip install "setuptools<70"
```
 
**Key pinned versions**:
- `pyscenic==0.12.1`
- `numpy==1.23.5`
- `scipy==1.10.1`
- `numba==0.57.1`
- Use `anndata` not `scanpy` in scripts to avoid matplotlib C++ conflicts
**Smoke test**:
```bash
module load gcc-11.2.0-gcc-8.5.0
/scratch/easmit31/conda_envs/pyscenic_final/bin/python -c "
import pandas as pd, numpy as np
from pyscenic.aucell import aucell
from pyscenic.utils import GeneSignature
data = pd.DataFrame(np.random.poisson(2, (200, 500)).astype(float),
                    columns=[f'gene{i}' for i in range(500)])
sig = GeneSignature(name='test', gene2weight={f'gene{i}': 1.0 for i in range(50)})
scores = aucell(data, [sig], num_workers=1)
print('✅ PASS' if scores['test'].mean() > 0 else '❌ FAIL')
"
```
 
---
 
## Design Notes
 
**Why liftover before GRN?** No macaque cisTarget databases exist. The aertslab only provides human (hg38), mouse, and fly databases. Liftover to human gene space before GRNBoost is the standard approach for NHP data.
 
**Why Python API for AUCell?** The `pyscenic aucell` CLI has a bug where it calls `adata.X.todense()` on HDF5 Dataset objects, causing a crash. The Python API with explicit dense conversion is the reliable alternative.
 
**Why mixed-effects models?** Cells from the same animal are not independent. Mixed-effects models with `animal_id` as random intercept correctly account for pseudoreplication.
 
