import argparse
import scanpy as sc
import scipy.sparse as sp

ap = argparse.ArgumentParser(description="Fix adult-only h5ad so .X is csr_matrix")
ap.add_argument("--prefix", required=True, help="Dataset prefix (e.g. gaba, astros)")
args = ap.parse_args()

in_file  = f"/scratch/easmit31/GRN_copy/scenic/h5ad_files/{args.prefix}_adults_allLv.h5ad"
out_file = f"/scratch/easmit31/GRN_copy/scenic/h5ad_files/{args.prefix}_adults_allLv_fixed.h5ad"

print(f"[fix_adults] Reading: {in_file}")
adata = sc.read_h5ad(in_file)   # full load into memory

# force into CSR sparse if needed
if not sp.issparse(adata.X):
    from scipy.sparse import csr_matrix
    adata.X = csr_matrix(adata.X)

print(f"[fix_adults] Cells={adata.n_obs}, Genes={adata.n_vars}, Matrix={type(adata.X)}")
adata.write(out_file)
print(f"[fix_adults] Wrote fixed file: {out_file}")
