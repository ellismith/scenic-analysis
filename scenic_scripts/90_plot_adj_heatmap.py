#!/usr/bin/env python
import os, argparse, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import loompy

# -----------------------------
# Parse command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Plot adjacency heatmap from SCENIC adj.csv and loom")
parser.add_argument("--adj", required=True, help="Adjacency CSV file from GRNBoost/SCENIC")
parser.add_argument("--loom", required=True, help="Loom file with gene order")
parser.add_argument("--tf-list", required=True, help="Text file with TFs, one per line")
parser.add_argument("--outd", required=True, help="Output directory")
args = parser.parse_args()
os.makedirs(args.outd, exist_ok=True)

# -----------------------------
# 1) Load gene order from loom
# -----------------------------
with loompy.connect(args.loom, "r") as ds:
    genes = list(map(str, ds.ra["Gene"]))
gene_index = {g: i for i, g in enumerate(genes)}
print(f"[info] genes in loom: {len(genes)}")

# -----------------------------
# 2) Load TF list and overlap
# -----------------------------
with open(args.tf_list) as f:
    tfs_all = [ln.strip() for ln in f if ln.strip()]
tfs = [t for t in tfs_all if t in gene_index]
tf_index = {t: i for i, t in enumerate(tfs)}
print(f"[info] TFs overlapping genes: {len(tfs)}")

# -----------------------------
# 3) Read adj.csv and map
# -----------------------------
adj = pd.read_csv(args.adj)
cols = {c.lower(): c for c in adj.columns}
TFcol = cols.get("tf") or cols.get("regulator") or list(adj.columns)[0]
TGcol = cols.get("target") or cols.get("gene") or list(adj.columns)[1]
IMcol = cols.get("importance") or cols.get("weight") or list(adj.columns)[2]
print(f"[info] using columns: TF={TFcol} target={TGcol} importance={IMcol}")

# filter to known TFs and genes
adj = adj[adj[TFcol].isin(tf_index) & adj[TGcol].isin(gene_index)]
print(f"[info] edges after filtering: {len(adj)}")
if adj.empty:
    raise SystemExit("[error] no edges after filtering; check adj.csv headers or TF list")

# map to integer indices
r = adj[TFcol].map(tf_index).to_numpy()
c = adj[TGcol].map(gene_index).to_numpy()
d = adj[IMcol].astype(np.float32).to_numpy()

# -----------------------------
# 4) Build sparse & dense matrix
# -----------------------------
nR, nC = len(tfs), len(genes)
M = coo_matrix((d, (r, c)), shape=(nR, nC)).tocsr().astype(np.float32)
D = M.toarray()
np.save(os.path.join(args.outd, "adjacency_dense.npy"), D)

# -----------------------------
# 5) Plot heatmap
# -----------------------------
vmax = np.percentile(D[D > 0], 99.5) if (D > 0).any() else 1.0
fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(np.log1p(D), aspect="auto")
ax.set_title("TF × Gene adjacency (importance, log1p)")
ax.set_xlabel(f"Genes (n={nC})")
ax.set_ylabel(f"TFs (n={nR})")
cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("log1p(importance)")
plt.tight_layout()
fig.savefig(os.path.join(args.outd, "adjacency_heatmap_full.png"), dpi=200)
plt.close(fig)

# -----------------------------
# 6) Row/col sums
# -----------------------------
row_sum = np.asarray(M.sum(axis=1)).ravel()
col_sum = np.asarray(M.sum(axis=0)).ravel()
pd.DataFrame({"TF": tfs, "row_sum_importance": row_sum}).to_csv(
    os.path.join(args.outd, "tf_row_sums.tsv"), sep="\t", index=False
)
pd.DataFrame({"Gene": genes, "col_sum_importance": col_sum}).to_csv(
    os.path.join(args.outd, "gene_col_sums.tsv"), sep="\t", index=False
)

# -----------------------------
# 7) Histograms
# -----------------------------
for arr, name in [(row_sum, "tf_row_sums"), (col_sum, "gene_col_sums"), (d, "edge_importances")]:
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(np.log1p(arr), bins=100)
    ax.set_title(f"{name} (log1p)")
    plt.tight_layout()
    fig.savefig(os.path.join(args.outd, f"{name}_hist.png"), dpi=150)
    plt.close(fig)

# -----------------------------
# 8) Save orders
# -----------------------------
with open(os.path.join(args.outd, "tf_order.txt"), "w") as f:
    f.write("\n".join(tfs))
with open(os.path.join(args.outd, "gene_order.txt"), "w") as f:
    f.write("\n".join(genes))

print("[ok] wrote:", os.path.join(args.outd, "adjacency_heatmap_full.png"))
