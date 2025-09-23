import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_auc_heatmap(adata, groupby="louvain", max_regulons=50, out_png="heatmap.png"):
    """
    Plot AUCell regulon activity heatmap grouped by cluster.

    Parameters
    ----------
    adata : AnnData
        Must have adata.obsm["AUCell"] with regulon activity
    groupby : str
        obs column to group by (default: 'louvain')
    max_regulons : int
        Number of top variable regulons to keep
    out_png : str
        Output filename
    """

    if "AUCell" not in adata.obsm:
        raise ValueError("adata.obsm['AUCell'] is missing")

    mat = adata.obsm["AUCell"]

    # --- debug info ---
    print(f"[plot] Input matrix shape: {mat.shape}")
    print(f"[plot] Any NaN? {mat.isna().any().any()}")
    print(f"[plot] Min={mat.min().min()}, Max={mat.max().max()}")

    # select top variable regulons
    variances = mat.var(axis=0)
    top = variances.sort_values(ascending=False).head(max_regulons).index
    sub = mat[top]

    # aggregate by group
    if groupby not in adata.obs:
        raise ValueError(f"Groupby column '{groupby}' not found in obs")
    grouped = sub.groupby(adata.obs[groupby]).mean()

    # --- more debug info ---
    print(f"[plot] After subsetting: {sub.shape}")
    print(f"[plot] After grouping: {grouped.shape}")
    print(f"[plot] Grouped NaN count: {grouped.isna().sum().sum()}")

    # if all NaN, bail out with message
    if grouped.isna().all().all():
        print("[plot] ❌ All NaN values after grouping — nothing to plot")
        return

    plt.figure(figsize=(0.25*grouped.shape[1]+4, 0.25*grouped.shape[0]+4))
    sns.heatmap(
        grouped,
        cmap="viridis",
        cbar_kws={"label": "AUCell score"},
        linewidths=0.5,
        linecolor="gray"
    )
    plt.title("AUCell regulon activity grouped by " + groupby)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"[plot] Saved heatmap to {out_png} (shape {grouped.shape})")
