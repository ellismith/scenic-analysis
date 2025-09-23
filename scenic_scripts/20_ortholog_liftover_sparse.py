#!/usr/bin/env python3
"""
20_ortholog_liftover_sparse.py
Map macaque genes to one-to-one human orthologs and collapse counts to
human gene symbols, producing a liftover AnnData (.h5ad).
This version includes progress logging during the collapse step.
"""

import sys, time
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

if len(sys.argv) != 4:
    sys.exit("Usage: python 20_ortholog_liftover_sparse.py <input.h5ad> <mapping.csv> <output.h5ad>")

input_h5ad, mapping_csv, out_h5ad = sys.argv[1:]

print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Input h5ad: {input_h5ad}", flush=True)
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Mapping CSV: {mapping_csv}", flush=True)
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Output h5ad: {out_h5ad}", flush=True)

# Load data
adata = sc.read_h5ad(input_h5ad)
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Loaded AnnData: {adata.shape}", flush=True)

# Load mapping
map_df = pd.read_csv(mapping_csv)
map_df = map_df[map_df["Human homology type"] == "ortholog_one2one"]
map_dict = dict(zip(map_df["Gene stable ID"], map_df["Human gene name"]))

adata.var["human_gene_name"] = adata.var_names.map(map_dict)
adata = adata[:, adata.var["human_gene_name"].notna()]
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Kept {adata.n_vars} genes with 1:1 orthologs", flush=True)

# Collapse to human genes
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Collapsing counts to human gene names (sparse)...", flush=True)
t0 = time.time()
mat = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

collapsed = {}
groups = adata.var.groupby("human_gene_name").indices
for i, (hgene, idxs) in enumerate(groups.items()):
    collapsed[hgene] = sp.csr_matrix(mat[:, list(idxs)].sum(axis=1))
    if i % 500 == 0:
        print(f"  Processed {i}/{len(groups)} genes...", flush=True)

human_mat = sp.hstack([collapsed[g] for g in collapsed.keys()]).tocsr()
adata_out = sc.AnnData(X=human_mat, obs=adata.obs.copy(), var=pd.DataFrame(index=list(collapsed.keys())))
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Collapse finished in {time.time()-t0:.1f} sec. Final shape: {adata_out.shape}", flush=True)

# Save
adata_out.write(out_h5ad)
print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Wrote {out_h5ad}", flush=True)
