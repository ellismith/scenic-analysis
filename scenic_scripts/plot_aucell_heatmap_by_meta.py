#!/usr/bin/env python3
import argparse, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def read_meta(path, meta_key, meta_cell_col):
    # Try index_col=0 (Scanpy export) then fallback
    try:
        meta = pd.read_csv(path, index_col=0)
        if meta_key not in meta.columns:
            raise ValueError
        meta.index.name = meta_cell_col
        return meta[[meta_key]]
    except Exception:
        meta = pd.read_csv(path)
        if meta_cell_col not in meta.columns or meta_key not in meta.columns:
            raise SystemExit(f"Metadata must have columns: {meta_cell_col} and {meta_key}")
        meta = meta.set_index(meta_cell_col)
        return meta[[meta_key]]

def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]

def main():
    ap = argparse.ArgumentParser(description="AUCell heatmap ordered by metadata (no clustering), with separators and sidebar.")
    ap.add_argument("--input", required=True, help="auc_mtx.csv")
    ap.add_argument("--meta", required=True, help="metadata CSV (index=cell, has meta-key)")
    ap.add_argument("--meta-key", required=True, help="column in metadata to order by (e.g. louvain)")
    ap.add_argument("--meta-cell-col", default="cell", help="name for the cell id (if not index)")
    ap.add_argument("--out", required=True, help="output PNG")
    ap.add_argument("--n-regs", type=int, default=30, help="number of regulons (most variable) to plot")
    ap.add_argument("--sample-per-group", type=int, default=30, help="balanced sample per group")
    ap.add_argument("--group-order", default="", help="comma-separated explicit order (e.g. 5,0,17,...)")
    ap.add_argument("--seed", type=int, default=0, help="random seed")
    # formatting
    ap.add_argument("--sep-color", default="white", help="separator line color")
    ap.add_argument("--sep-width", type=float, default=1.6, help="separator line width")
    ap.add_argument("--label-fontsize", type=int, default=12, help="fontsize for group labels")
    ap.add_argument("--reg-fontsize", type=int, default=9, help="fontsize for regulon labels")
    ap.add_argument("--title", default="", help="custom plot title")
    ap.add_argument("--sidebar", action="store_true", help="draw a left color strip indicating groups")
    ap.add_argument("--sidebar-legend", action="store_true", help="add legend for the sidebar colors")
    args = ap.parse_args()

    # Read AUCell (cells x regulons), first column = cell ID
    auc = pd.read_csv(args.input)
    cell_col = auc.columns[0]
    auc = auc.set_index(cell_col)

    # Read metadata
    meta = read_meta(args.meta, args.meta_key, args.meta_cell_col)

    # Intersect cells
    common = auc.index.intersection(meta.index)
    if len(common) == 0:
        raise SystemExit("No overlap between AUCell cells and metadata cells.")
    auc = auc.loc[common]
    meta = meta.loc[common]

    # Group order
    groups = meta[args.meta_key].astype(str)
    if args.group_order:
        wanted = [g.strip() for g in args.group_order.split(",") if g.strip()]
        order = [g for g in wanted if g in set(groups)]
        # append any unseen groups at the end
        order += [g for g in sorted(groups.unique(), key=natural_key) if g not in order]
    else:
        order = sorted(groups.unique(), key=natural_key)

    # Balanced sampling per group
    rng = np.random.default_rng(args.seed)
    idx_list = []
    for g in order:
        g_cells = meta.index[groups == g].tolist()
        take = min(args.sample_per_group, len(g_cells))
        if take > 0:
            pick = rng.choice(g_cells, size=take, replace=False)
            idx_list.extend(pick.tolist())

    if len(idx_list) == 0:
        raise SystemExit("No cells selected for plotting; try increasing --sample-per-group.")

    sub_auc = auc.loc[idx_list]
    sub_meta = meta.loc[idx_list]

    # Sort rows by group then by cell id
    row_order = sorted(sub_auc.index, key=lambda c: (order.index(str(sub_meta.loc[c, args.meta_key])), c))
    sub_auc = sub_auc.loc[row_order]
    sub_meta = sub_meta.loc[row_order]

    # Pick most variable regulons (within the subset)
    variances = sub_auc.var(axis=0)
    regs = list(variances.sort_values(ascending=False).head(args.n_regs).index)
    mat = sub_auc[regs].to_numpy(dtype=float)

    # Scale each regulon to [0,1] for contrast
    col_min = mat.min(axis=0, keepdims=True)
    col_max = mat.max(axis=0, keepdims=True)
    denom = np.where((col_max - col_min) == 0, 1.0, col_max - col_min)
    matn = (mat - col_min) / denom

    # Compute group sizes & boundaries
    group_labels = sub_meta[args.meta_key].astype(str).values
    sizes = [np.sum(group_labels == g) for g in order]
    bounds = np.cumsum([s for s in sizes if s > 0])

    # Figure size
    h, w = matn.shape
    fig_w = min(26, 2 + 0.36 * w)   # a bit wider labels
    fig_h = min(22, 2 + 0.12 * h)
    if args.sidebar:
        fig_w += 1.4  # room for the sidebar
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Axes layout: optional sidebar + main heatmap
    if args.sidebar:
        gs = fig.add_gridspec(ncols=2, nrows=1, width_ratios=[0.04, 1.0], wspace=0.05)
        ax_bar = fig.add_subplot(gs[0,0])
        ax = fig.add_subplot(gs[0,1])
    else:
        ax = fig.add_subplot(111)
        ax_bar = None

    # Main heatmap
    im = ax.imshow(matn, aspect="auto", interpolation="nearest")
    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("AUCell (scaled 0–1 per regulon)")

    # Regulon labels
    ax.set_xticks(np.arange(w))
    ax.set_xticklabels(regs, rotation=60, ha="right", fontsize=args.reg_fontsize)
    # Y ticks: group names at midpoints
    y_ticks = []
    y_labels = []
    start = 0
    for g, s in zip(order, sizes):
        if s > 0:
            y_ticks.append(start + s/2)
            y_labels.append(str(g))
            start += s
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=args.label_fontsize, fontweight="bold")

    # Separator lines
    for b in bounds[:-1]:
        ax.axhline(b - 0.5, color=args.sep_color, linewidth=args.sep_width)

    # Title
    title = args.title or f"AUCell heatmap ordered by {args.meta_key} — {h} cells × {w} regulons"
    ax.set_title(title, fontsize=max(12, args.label_fontsize+2), pad=10)

    # Optional sidebar color strip (acts like a bracket)
    if ax_bar is not None:
        # map each group to a color
        uniq = [g for g in order if g in set(group_labels)]
        # repeat colors if too many groups
        base_colors = plt.cm.tab20.colors
        colors = [base_colors[i % len(base_colors)] for i in range(len(uniq))]
        cmap = ListedColormap(colors)
        # build a column of group indices
        g_to_idx = {g:i for i,g in enumerate(uniq)}
        col = np.array([g_to_idx[g] for g in group_labels], dtype=int)[:, None]  # (rows, 1)
        ax_bar.imshow(col, aspect="auto", cmap=cmap, interpolation="nearest")
        ax_bar.set_xticks([]); ax_bar.set_yticks([])
        for b in bounds[:-1]:
            ax_bar.axhline(b - 0.5, color=args.sep_color, linewidth=args.sep_width)

        if args.sidebar_legend:
            # tiny legend to map color → group
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=colors[i], edgecolor='none', label=str(g)) for i,g in enumerate(uniq)]
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, title=args.meta_key)

    plt.tight_layout()
    plt.savefig(args.out, dpi=220)
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
