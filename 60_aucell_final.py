#!/usr/bin/env python3
"""
60_aucell_patched2.py
Run AUCell with proper gene ranking, write regulon scores (no AUC_ prefix).
"""
import sys
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from pyscenic.aucell import aucell
from pyscenic.utils import GeneSignature

if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} <expr.h5ad> <regulons.csv> <out.h5ad>")

expr_h5ad, regulons_csv, out_h5ad = sys.argv[1:]

print("👉 Loading AnnData:", expr_h5ad, flush=True)
adata = sc.read_h5ad(expr_h5ad, backed=None)

if not isinstance(adata.X, sp.csr_matrix):
    adata.X = sp.csr_matrix(adata.X)

print("👉 Building expression DataFrame", flush=True)
ex_mtx = pd.DataFrame.sparse.from_spmatrix(
    adata.X,
    index=adata.obs_names,
    columns=adata.var_names,
)

print("👉 Loading regulons:", regulons_csv, flush=True)
df = pd.read_csv(regulons_csv)
if not {"regulon", "target"}.issubset(df.columns):
    sys.exit("CSV must have columns: regulon,target")

# Build regulon signatures
regulon_sigs = []
for reg, group in df.groupby("regulon"):
    genes = list(group["target"].dropna().unique())
    regulon_sigs.append(GeneSignature(name=reg, gene2weight={g: 1.0 for g in genes}))

print(f"👉 Built {len(regulon_sigs)} regulon signatures", flush=True)

print("👉 Running AUCell (no pre-ranking)…", flush=True)
scores = aucell(ex_mtx, regulon_sigs, num_workers=8)

print("👉 Attaching results (no AUC_ prefix)", flush=True)
for sig in regulon_sigs:
    adata.obs[sig.name] = scores[sig.name]

print("👉 Writing output:", out_h5ad, flush=True)
adata.write(out_h5ad, compression="gzip")
print("✅ Done. Wrote:", out_h5ad, flush=True)
