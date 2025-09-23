import scanpy as sc
import pandas as pd

print("Loading big AnnData slice in backed mode...")
adata_backed = sc.read_h5ad(
    "/scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_adults_allLv_ready.h5ad",
    backed="r"
)

# Subset by indices without loading everything
obs_idx = list(range(1000))   # first 1000 cells
var_idx = list(range(100))    # first 100 genes
adata = adata_backed[obs_idx, var_idx].to_memory()
adata.write("test_small.h5ad")
print("Wrote test_small.h5ad")

map_df = pd.DataFrame({
    "Gene stable ID": adata.var_names,
    "Human gene name": [f"HG{i}" for i in range(len(adata.var_names))],
    "Human homology type": ["ortholog_one2one"]*len(adata.var_names),
})
map_df.to_csv("test_map.csv", index=False)
print("Wrote test_map.csv")
