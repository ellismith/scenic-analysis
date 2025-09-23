import sys
import scanpy as sc
from scipy import sparse

if len(sys.argv) != 3:
    sys.exit("Usage: python fix_h5ad_sparse.py <in.h5ad> <out.h5ad>")

in_path, out_path = sys.argv[1], sys.argv[2]
print(f"[fix_h5ad_sparse] Reading {in_path}")
adata = sc.read_h5ad(in_path, backed=None)

if not sparse.isspmatrix_csr(adata.X):
    print("[fix_h5ad_sparse] Converting .X to csr_matrix")
    adata.X = sparse.csr_matrix(adata.X)

adata.write(out_path)
print(f"[fix_h5ad_sparse] Wrote {out_path}")
