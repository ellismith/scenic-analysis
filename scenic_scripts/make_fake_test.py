import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

print("Creating fake test AnnData...")
X = sp.random(1000, 100, density=0.05, format="csr", dtype=np.float32)  # 1000 cells × 100 genes
obs = pd.DataFrame(index=[f"cell{i}" for i in range(1000)])
var = pd.DataFrame(index=[f"gene{i}" for i in range(100)])

adata = sc.AnnData(X=X, obs=obs, var=var)
adata.write("test_small.h5ad")
print("Wrote test_small.h5ad")

map_df = pd.DataFrame({
    "Gene stable ID": adata.var_names,
    "Human gene name": [f"HG{i//2}" for i in range(len(adata.var_names))],  # collapse 2 genes → 1 human
    "Human homology type": ["ortholog_one2one"]*len(adata.var_names),
})
map_df.to_csv("test_map.csv", index=False)
print("Wrote test_map.csv")
