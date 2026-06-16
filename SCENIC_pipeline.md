# pySCENIC GRN Pipeline — Complete Guide
## Macaque Brain Aging

Goal: Identify transcription factor (TF) regulons (TF and its coexpressed target genes), and in what brain regions / cell types that...
1. Show age- and sex-related changes in regulon activity 
2. Show high (or low) inter-individual variability
3. Show inter-individual variability increases, decreases, or no changes with age

How: run the pySCENIC gene regulatory network pipeline on macaque snRNA-seq data

References: 
SCENIC website (links to their papers): https://scenic.aertslab.org/
Documentation: https://github.com/aertslab/SCENICprotocol/tree/master/notebooks

Preprint of the dataset: https://pubmed.ncbi.nlm.nih.gov/41279774/

---

## Getting Started

### Source data location: contains the h5ad files with the RNAseq data

U01 dataset (55 animals, 11 brain regions)
```
/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/
```

### Scripts location
```
/scratch/easmit31/GRN/scenic/scenic_scripts/
https://github.com/ellismith/scenic-analysis/
```

### Reference databases 
```
/scratch/easmit31/data/cisTarget_databases/
/scratch/easmit31/GRN/scenic/allTFs_hg38.txt
```

---

## Files to Copy

Copy the following to your own scratch directory. Note: Some of these file paths are currently hardcoded in the scripts, so you'll need to either edit the path or we can just make it flexible.

### Reference databases (copy)
```bash
cp -r /scratch/easmit31/data/cisTarget_databases /scratch/<yourasuid>/data/

# TF list
cp /scratch/easmit31/GRN/scenic/allTFs_hg38.txt /scratch/<yourasuid>/GRN/scenic/
```

### Conda environment yaml
```bash
cp /scratch/easmit31/conda_envs/pyscenic_env.yaml /scratch/<yourasuid>/
```


### Create output directories
```bash
mkdir -p /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs
mkdir -p /scratch/<yourasuid>/GRN/scenic/results
mkdir -p /scratch/<yourasuid>/GRN/scenic/scenic_scripts/logs
```

---

## Environment Setup

### Build the conda environment
```bash
module load mamba/latest
mamba env create -p /scratch/<yourasuid>/conda_envs/pyscenic_final \
  -f /path/to/pyscenic_env.yaml

# Fix setuptools (required — pyscenic 0.12.1 needs pkg_resources)
/scratch/<yourasuid>/conda_envs/pyscenic_final/bin/pip install "setuptools<70"
```

### Always load gcc before running anything
```bash
module load gcc-11.2.0-gcc-8.5.0
```
Without this you get `GLIBCXX_3.4.30 not found` errors. It's already included in all `.sh` scripts.

### Pinned versions 
| Package | Version | Why pinned |
|---------|---------|------------|
| pyscenic | 0.12.1 | Stable API |
| numpy | 1.23.5 | GLIBCXX compatibility |
| scipy | 1.10.1 | GLIBCXX compatibility |
| numba | 0.57.1 | GLIBCXX compatibility |

### Verify environment works
```bash
module load gcc-11.2.0-gcc-8.5.0
/scratch/<yourasuid>/conda_envs/pyscenic_final/bin/python -c "
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

## Pipeline Steps


---

### Step 10 — Select cells 
**⚠️ NOTE: THIS STEP CAN BE COLLAPSED INTO THE NEXT ONE... we shouldn't need to re-save the h5ad files, so just filter the cells by to exclude the 2 infants (to do this I just use min age of 1 year to only include the cells from the adults)**

**Script**: `10_select_cells.sh`
**Purpose**: Filter source h5ad to a single cell type, adults only (age ≥ 1 year) -- just collapse into next step. 
**Mem**: 64G (this is a guess, edit this as you go through)

```bash
cd /scratch/<yourasuid>/GRN/scenic/scenic_scripts
sbatch 10_select_cells.sh [args]
```

**Output**: `h5ad_adult_files/adult_grn_<celltype>_selected.h5ad`

**Output structure**:
```
shape: for microglia, ~158k cells × 29,060 macaque genes (varies by cell type, there are 1.4M glutamatergic neurons)
obs: age, animal_id, sex, louvain, ct_louvain, region, brain_id, ...
var: _index (Ensembl ID), ensembl_gene_id, external_gene_name, chr, bp1, bp2, mt, strand
X: raw integer counts
```

Other note, OPCs and oligodendrocytes are combined in the h5ad file, but the OPCs are louvains 12 and 13, the rest are oligodendrocytes.

**Check**: `var['chr']` must be present — it's used for sex-chr filtering in step 74.
We may want to move the filtering to remove genes on sex chromosomes upstream though. 
29,060 total, 1,344 sex-chr, 27,716 autosomal

---

### Step 20 — Ortholog liftover
**Script**: `20_ortholog_liftover.sh` (calls `20_ortholog_liftover_sparse.py`)
**Purpose**: Map macaque Ensembl IDs → human HGNC symbols (1:1 orthologs only)
**Mem**: 512G (this is a guess, edit this as you go through)

```bash
sbatch 20_ortholog_liftover.sh \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_selected.h5ad \
  /path/to/ensembl_macaque_human_mapping.csv \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_liftover.h5ad
```

**⚠️ NOTE: Edit to use all cells instead of subsetting? **: Edit `20_ortholog_liftover_sparse.py`:
```python
# Find this argument:
parser.add_argument("--downsample", type=int, default=100000, ...)

# Change default to 0, and update the condition from:
if args.downsample and adata.n_obs > args.downsample:
# To:
if args.downsample and args.downsample > 0 and adata.n_obs > args.downsample:
```

**Output**: `h5ad_adult_files/adult_grn_<celltype>_liftover.h5ad`

**Output structure**:
```
shape: up to N cells × ~12,788 human genes
obs: same as selected (all metadata preserved)
var: human_gene_name (index), n_mm_geneid
X: raw integer counts (collapsed from macaque → human)
```

I think the right approach is liftover before GRNBoost, but let me know if this seems wrong. 

We can also run the pipeline on a subset of highly variable genes rather than all genes after liftover.

---

### Step 30 — Convert to loom
**Script**: `30_convert_to_loom.sh`
**Purpose**: Convert h5ad to loom format required by GRNBoost2
**Mem**: 256G (this is a guess, edit this as you go through)

```bash
sbatch 30_convert_to_loom.sh \
  --adata_in /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_liftover.h5ad \
  --loom_out /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>.loom
```

**Output**: `adult_grn_outputs/adult_grn_<celltype>.loom`

**Output structure**:
```
shape: N cells × 12,788 genes
row attrs: Gene (human HGNC symbols)
col attrs: CellID, nGene, nUMI
matrix: float32, stored as genes × cells (transposed from h5ad)
```

**Check**:
```bash
python3 -c "
import loompy
ds = loompy.connect('adult_grn_<celltype>.loom', 'r')
print(ds.shape, list(ds.ra.keys()), list(ds.ca.keys()))
ds.close()
"
```

---

### Step 40 — GRN inference (GRNBoost2)
**Script**: `40_run_grn.sh`
**Purpose**: Infer TF-gene co-expression network
**Partition**: `htc` | **Mem**: 200G | **CPUs**: 20

```bash
sbatch 40_run_grn.sh \
  --loom_in /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>.loom \
  --tfs /scratch/<yourasuid>/GRN/scenic/allTFs_hg38.txt \
  --out /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>_adj.csv
```

**Output**: `adult_grn_outputs/adult_grn_<celltype>_adj.csv`

**Output structure**:
```
columns: TF, target, importance
rows: ~1M TF-gene pairs
importance range: 0 – ~15
```

**Check**:
```bash
head -5 adult_grn_<celltype>_adj.csv
wc -l adult_grn_<celltype>_adj.csv  # expect ~1M rows
```

---

### Step 50 — cisTarget pruning
**Script**: `50_run_ctx.sh`
**Purpose**: Prune adjacencies using motif databases → regulons
**Mem**: 64G | **CPUs**: 8  (this is a guess, edit this as you go through)

```bash
sbatch 50_run_ctx.sh \
  --adj /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>_adj.csv \
  --expression_mtx_fname /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>.loom \
  --db_glob "/scratch/<yourasuid>/data/cisTarget_databases/mc_v10_clust_ranking/human/hg38_*.rankings.feather" \
  --annotations_fname /scratch/<yourasuid>/data/cisTarget_databases/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --output /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>_regulons.csv \
  --workers 8
```

**Output**: `adult_grn_outputs/adult_grn_<celltype>_regulons.csv`

**Output structure**: Raw cisTarget CSV with 3-line header. Contains TF, MotifID, AUC, NES, TargetGenes (as tuple list). Parsed by step 56.

**Check**:
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('adult_grn_<celltype>_regulons.csv', skiprows=3, header=None)
print(f'{len(df)} regulons')
print(df[0].value_counts().head(10))
"
```
Expect: 1,000–3,000 regulons. Top TFs should be biologically sensible for your cell type.

---

### Step 56 — Convert regulons
**Script**: `56_convert_regulons.py`
**Purpose**: Parse raw cisTarget CSV → simple `regulon,target` pairs
**Runtime**: ~1 min (run interactively or in a short job)

```bash
/scratch/<yourasuid>/conda_envs/pyscenic_final/bin/python \
  /scratch/<yourasuid>/GRN/scenic/scenic_scripts/56_convert_regulons.py \
  /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>_regulons.csv \
  /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>_regulons_simple.csv
```

**Output**: `adult_grn_outputs/adult_grn_<celltype>_regulons_simple.csv`

**Output structure**:
```
columns: regulon, target
rows: ~64k regulon-target pairs (varies)
unique regulons: ~297 (after collapsing multiple motifs per TF)
mean targets per regulon: ~215
```

**Check**:
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('adult_grn_<celltype>_regulons_simple.csv')
sizes = df.groupby('regulon').size().sort_values(ascending=False)
print(f'{df.regulon.nunique()} unique regulons')
print(sizes.head(10))  # check for artifact regulons with >5000 targets
"
```

---

### Step 60 — AUCell scoring
AUCell (Area Under the Curve for Cell scoring) 

For each cell, it ranks all genes from most to least expressed. Then for each regulon (a TF + its target genes), it asks: are the regulon's target genes enriched near the top of that ranking? 

It computes an AUC — the area under the recovery curve — which measures how much the target genes cluster at the top of the expression ranking compared to random.

AUC score is a number between 0 and 1. 0 means the regulon's target genes are no more highly expressed in that cell than you'd expect by chance. Near 1 means essentially all of the regulon's targets are among the most highly expressed genes in that cell. 

A high AUC score for a given TF regulon in a cell means that TF's transcriptional program is active in that cell -- whether the downstream gene program the TF controls is turned on. So, a TF can be expressed but not active, or active at low expression levels.

Why not just use TF expression directly? TF mRNA levels are noisy and don't always reflect activity. AUCell captures the coordinated expression of the whole regulon, which is a more robust signal of actual TF activity.

**Script**: `60_aucell_final.sh`
**Purpose**: Score each cell for regulon activity
**Partition**: `htc` (≤200k cells) or `highmem` (>200k) | **Mem**: 200G  (this is a guess, edit this as you go through)

```bash
sbatch 60_aucell_final.sh \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_liftover.h5ad \
  /scratch/<yourasuid>/GRN/scenic/adult_grn_outputs/adult_grn_<celltype>_regulons_simple.csv \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_aucell.h5ad
```

**For large cell types (>200k cells)**: Before, I tried this chunking approach, but if there's a way to avoid that, it would be better to do all at once:
```bash
# Chunk into 50k pieces
python3 chunk_h5ad.py \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_liftover.h5ad \
  50000 \
  adult_grn_<celltype>

# Submit AUCell per chunk using chunk_job.sh
# Then merge:
bash merge_aucell_chunks.sh <celltype>
```

**Output**: `h5ad_adult_files/adult_grn_<celltype>_aucell.h5ad`

**Output structure**:
```
shape: N cells × same vars as liftover
obs: all original metadata + one column per regulon (bare TF name, NO prefix)
     e.g. STAT1, IRF1, PPARG, ... (not AUC_STAT1)
var: same as liftover (human_gene_name, n_mm_geneid)
X: same as liftover (unchanged)
score range: 0.0 – ~0.5
nonzero per regulon: ~15–60% of cells
```

**⚠️ Critical check before proceeding**:
```bash
module load gcc-11.2.0-gcc-8.5.0
/scratch/<yourasuid>/conda_envs/pyscenic_final/bin/python -c "
import h5py, numpy as np
with h5py.File('adult_grn_<celltype>_aucell.h5ad', 'r') as f:
    meta = {'age','sex','animal_id','louvain','ct_louvain','region','_index',
            'brain_id','hemisphere','n_umi','total_counts','total_counts_mt'}
    cols = [k for k in f['obs'].keys()
            if k not in meta and isinstance(f['obs'][k], h5py.Dataset)]
    print(f'{len(cols)} regulons found')
    scores = f['obs'][cols[0]][:]
    print(f'Score range: {scores.min():.4f} - {scores.max():.4f}')
    print(f'Non-zero: {np.sum(scores>0)}/{len(scores)}')
    if scores.max() == 0:
        print('❌ ALL ZEROS — do not proceed, debug AUCell')
    else:
        print('✅ Scores look valid')
"
```
If all scores are zero, stop. See [Known Issues](#known-issues-and-fixes).

---

### Step 70 — Global mixed-effects models
Note: This might not be necessary..
**Script**: `70_regulon_mixed_models.sh`
**Purpose**: Model regulon activity ~ age + sex across all cells
**Partition**: `htc` |  **Mem**: 64G

```bash
sbatch 70_regulon_mixed_models.sh \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_aucell.h5ad \
  /scratch/<yourasuid>/GRN/scenic/results/adult_grn_<celltype>_mixed_models.tsv
```

**Output**: `results/adult_grn_<celltype>_mixed_models.tsv`

**Output structure**:
```
columns: regulon, term, coef, se, pval, random_effect_var, residual_var, icc, converged, qval
terms: Intercept, age_std, sex[T.M]
FDR correction: BH per term (excluding Intercept)
```

**Model**: `activity_std ~ age_std + sex` with random intercept per `animal_id`
Both age and activity are z-scored before modeling.

**Check**:
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('adult_grn_<celltype>_mixed_models.tsv', sep='\t')
age = df[df['term']=='age_std'].sort_values('qval')
print(f'Significant age regulons (q<0.05): {(age[\"qval\"]<0.05).sum()} / {len(age)}')
print(age[['regulon','coef','qval']].head(10).to_string(index=False))
"
```

---

### Step 74 — Cluster-level mixed-effects models
**Script**: `74_cluster_mixed_models.sh`
**Purpose**: Same models per louvain cluster, with sex-chr TF filtering
**Partition**: `htc` | **Time**: 2–4 hrs | **Mem**: 128G

```bash
sbatch 74_cluster_mixed_models.sh \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_aucell.h5ad \
  /scratch/<yourasuid>/GRN/scenic/h5ad_adult_files/adult_grn_<celltype>_selected.h5ad \
  /scratch/<yourasuid>/GRN/scenic/results/adult_grn_<celltype>_cluster_mixed_models.tsv
```

The second argument (selected h5ad) provides `var['chr']` for sex-chr filtering.

**Output**: `results/adult_grn_<celltype>_cluster_mixed_models.tsv`

**Output structure**:
```
columns: cluster, regulon, term, coef, se, pval, icc, n_cells, n_animals, qval
FDR correction: BH per term, across ALL clusters jointly
```

---

## Data Structure Reference

### Complete file map
```
/scratch/<yourasuid>/GRN/scenic/
├── allTFs_hg38.txt                         # 1,892 human TFs
├── h5ad_adult_files/
│   ├── adult_grn_<ct>_selected.h5ad        # Step 10 output
│   ├── adult_grn_<ct>_liftover.h5ad        # Step 20 output
│   └── adult_grn_<ct>_aucell.h5ad          # Step 60 output
├── adult_grn_outputs/
│   ├── adult_grn_<ct>.loom                 # Step 30 output
│   ├── adult_grn_<ct>_adj.csv              # Step 40 output
│   ├── adult_grn_<ct>_regulons.csv         # Step 50 output
│   └── adult_grn_<ct>_regulons_simple.csv  # Step 56 output
└── results/
    ├── adult_grn_<ct>_mixed_models.tsv              # Step 70 output
    └── adult_grn_<ct>_cluster_mixed_models.tsv      # Step 74 output

/scratch/<yourasuid>/data/cisTarget_databases/
├── mc_v10_clust_ranking/human/
│   ├── hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
│   └── hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
└── motif2tf/
    └── motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

### Data shape at each step (microglia reference)
| Step | File | Cells | Genes/Vars | Notes |
|------|------|-------|------------|-------|
| 10 | selected.h5ad | 158,201 | 29,060 | macaque genes |
| 20 | liftover.h5ad | 100,000* | 12,788 | human genes, *downsampled |
| 30 | .loom | 100,000 | 12,788 | float32 |
| 40 | _adj.csv | 1,062,195 rows | 3 cols | TF, target, importance |
| 50 | _regulons.csv | 2,051 regulons | — | raw cisTarget format |
| 56 | _regulons_simple.csv | 64,011 rows | 2 cols | regulon, target |
| 60 | _aucell.h5ad | 100,000 | 335 regulon scores in obs | |
| 70 | _mixed_models.tsv | 948 rows | 10 cols | 316 regulons × 3 terms |


---

## Known Issues and Fixes

### AUCell returns all-zero scores
**Cause**: Passing a sparse-backed pandas DataFrame to `aucell()`.
**Fix**: Already handled in `60_aucell_final.py` — it converts to dense with `.toarray()`. If you still get zeros, check that the expression matrix has non-zero values:
```python
print(adata.X[:5, :5].toarray())  # should NOT be all zeros
```

### `pkg_resources` ModuleNotFoundError
**Cause**: setuptools ≥70 removed `pkg_resources` from default namespace.
**Fix**: `pip install "setuptools<70"`

### `pyscenic ctx` picks up wrong binary
**Cause**: A conflicting pyscenic installation in `~/.local/`.
**Fix**: Use full path in `50_run_ctx.py`:
```python
cmd = ["/scratch/<yourasuid>/conda_envs/pyscenic_final/bin/pyscenic", "ctx", ...]
```

### GLIBCXX errors from numba/llvmlite
**Cause**: gcc module not loaded.
**Fix**: `module load gcc-11.2.0-gcc-8.5.0` before any pySCENIC work.

### Mixed models find no regulon columns
**Cause**: Scripts look for `AUC_` prefix but `60_aucell_final.py` writes bare TF names.
**Fix**: Already fixed in current scripts — they filter by `h5py.Dataset` type and exclude known metadata keys.

### `general` partition error on Sol
**Cause**: Sol renamed `general` to `public` in May 2026.
**Fix**: Use `--partition=public` in any script that had `--partition=general`.
