# SCENIC Gene Regulatory Network Analysis Pipeline

A comprehensive pipeline for inferring and analyzing gene regulatory networks (GRNs) using SCENIC/pySCENIC on single-cell RNA-seq data from macaque brain tissue.

## Overview

This pipeline processes single-cell expression data to:
- Infer transcription factor regulatory networks using GRNBoost2
- Calculate regulon activity scores with AUCell
- Perform comparative analysis across cell types and developmental stages

## Directory Structure

```
scenic/
├── h5ad_files/          # Expression data and regulon definitions
├── results/             # Statistical analysis outputs
├── figures_all/         # Full dataset visualizations
└── figures_subset/      # Subset dataset visualizations
```

## Main Data Files

### Expression Matrices
- `gaba_adults_allLv_pyscenic_output_merged.h5ad` (523K cells)
- `astros_adults_allLv_pyscenic_output_merged.h5ad` (416K cells)

### Regulon Definitions
- `gaba_allLv_regulons.csv` (331 regulons)
- `astros_allLv_regulons.csv` (214 regulons)

Format: `regulon,target` pairs defining TF-gene relationships

## Pipeline Workflow

### 1. Cell Selection
```bash
sbatch 10_select_cells.sh <input.h5ad> <output.h5ad> --ages adults,infants
```
Filter cells by age group (adults/infants), cell subclusters (select specific Louvain clusters to include), and brain region 

### 2. Ortholog Liftover
```bash
sbatch 20_ortholog_liftover.sh <input.h5ad> <output.h5ad> <ortholog_table.csv>
```
Convert macaque genes to human HGNC symbols

### 3. Convert to Loom
```bash
sbatch 30_convert_to_loom.sh <input.h5ad> <output.loom>
```
Prepare data for GRNBoost2

### 4. Network Inference (GRNBoost2)
```bash
sbatch 40_run_grnboost.sh <input.loom> <output_adj.csv>
```
Infer TF-target regulatory relationships

### 5. Motif Pruning (cisTarget)
```bash
sbatch 50_run_ctx.sh <input_adj.csv> <motif_db> <regulons.csv>
```
Refine regulons using motif enrichment

### 6. AUCell Scoring
Calculate regulon activity scores per cell using the generated regulon definitions

## Requirements

- Python 3.8+
- pySCENIC
- scanpy
- pandas, numpy
- SLURM cluster environment

## Key Results

- **Mixed-effects models**: Regulon activity variation across conditions
- **Variance analysis**: Cell type-specific regulatory differences
- **Visualizations**: UMAP projections, heatmaps, and comparative plots

## Usage Notes

- Use subset datasets for testing/development
- Full datasets require high-memory nodes (>180GB for GABAergic neurons)
- Regulon CSVs are the key upstream outputs for all downstream analyses

## Contact

For questions about this pipeline, please open an issue or contact elliot.a.smith@asu.edu.
