import sys
import scanpy as sc

if len(sys.argv) < 2:
    print("Usage: python summarize_h5ad.py <your_file.h5ad>")
    sys.exit(1)

f = sys.argv[1]
adata = sc.read_h5ad(f)

print("\n=== AnnData summary ===")
print(adata)

print("\n=== .obs (cell metadata) ===")
print(f"Shape: {adata.obs.shape}")
print("Columns:", list(adata.obs.columns)[:15], "..." if len(adata.obs.columns) > 15 else "")

print("\n=== .var (gene metadata) ===")
print(f"Shape: {adata.var.shape}")
print("Columns:", list(adata.var.columns)[:15], "..." if len(adata.var.columns) > 15 else "")

print("\n=== .uns (unstructured annotations) ===")
print("Keys:", list(adata.uns.keys()))

print("\n=== .obsm (cell-level matrices, e.g. PCA/UMAP) ===")
print("Keys:", list(adata.obsm.keys()))

print("\n=== .varm (gene-level matrices, e.g. loadings) ===")
print("Keys:", list(adata.varm.keys()))

print("\n=== .obsp (cell × cell pairwise matrices) ===")
print("Keys:", list(adata.obsp.keys()))

print("\n=== .varp (gene × gene pairwise matrices) ===")
print("Keys:", list(adata.varp.keys()))

print("\n=== .layers (alternative expression matrices) ===")
print("Keys:", list(adata.layers.keys()))

# Optional: peek into uns contents
for k in adata.uns.keys():
    print(f"  uns['{k}']: type={type(adata.uns[k])}")
