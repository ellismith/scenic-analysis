#!/usr/bin/env python3
"""
plot_aucell_heatmap.py — Flexible AUCell heatmap

Aggregation levels (--level):
  cell            raw single-cell scores (default)
  animal          pseudobulk mean per animal
  cluster         mean per louvain cluster
  animal_region   mean per animal x region

Regulon selection (--select-regs):
  variable_between   variance of animal means (inter-individual, default)
  variable_within    mean of per-animal variance (intra-individual)
  variable_cluster   variance across cluster means
  mean_activity      top mean AUCell score
  variable_cell      variance across all cells

Sort rows (--sort-by): any metadata column (age, region, sex, louvain, animal_id)
                       or 'mean_activity'

Examples:
  # Animal level, sorted by age, top between-individual variable regulons
  python plot_aucell_heatmap.py --input scores.csv --meta aucell.h5ad \\
    --level animal --sort-by age --ytick-cols animal_id,age \\
    --select-regs variable_between --n-regs 30 --cluster-regs --out out.png

  # Cell level, sorted by region, colored by region
  python plot_aucell_heatmap.py --input scores.csv --meta aucell.h5ad \\
    --level cell --sort-by region --color-by region \\
    --select-regs variable_between --n-regs 30 --n-cells 2000 --out out.png

  # Specific regulons
  python plot_aucell_heatmap.py --input scores.csv --meta aucell.h5ad \\
    --level animal --sort-by age --ytick-cols animal_id,age \\
    --regulons PPARG,IRF3,MAF,MITF --out out.png
"""
import argparse
import sys
import h5py
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

META_COLS = ['louvain', 'ct_louvain', 'age', 'sex', 'region', 'animal_id', 'brain_id']

def load_meta(h5ad_path):
    with h5py.File(h5ad_path, 'r') as f:
        cell_ids = f['obs']['_index'][:].astype(str)
        meta = pd.DataFrame(index=cell_ids)
        for col in META_COLS:
            if col not in f['obs']:
                continue
            grp = f['obs'][col]
            if isinstance(grp, h5py.Dataset):
                meta[col] = grp[:]
            elif 'codes' in grp and 'categories' in grp:
                cats = grp['categories'][:].astype(str)
                codes = grp['codes'][:]
                meta[col] = [cats[c] for c in codes]
    return meta

def select_regulons(df, meta, method, n, specified):
    if specified.strip():
        regs = [r.strip() for r in specified.split(',') if r.strip() in df.columns]
        if not regs:
            sys.exit("None of the specified regulons found.")
        return regs

    if method == 'variable_cell':
        return list(df.var().sort_values(ascending=False).head(n).index)

    if method == 'mean_activity':
        return list(df.mean().sort_values(ascending=False).head(n).index)

    if 'animal_id' not in meta.columns:
        print("Warning: animal_id not in metadata, falling back to variable_cell")
        return list(df.var().sort_values(ascending=False).head(n).index)

    animal_means = df.join(meta['animal_id']).groupby('animal_id')[df.columns].mean()

    if method == 'variable_between':
        return list(animal_means.var().sort_values(ascending=False).head(n).index)

    if method == 'variable_within':
        within = df.join(meta['animal_id']).groupby('animal_id')[df.columns].var().mean()
        return list(within.sort_values(ascending=False).head(n).index)

    if method == 'variable_cluster':
        grp_col = 'louvain' if 'louvain' in meta.columns else None
        if not grp_col:
            sys.exit("louvain not in metadata, cannot use variable_cluster")
        cluster_means = df.join(meta[grp_col]).groupby(grp_col)[df.columns].mean()
        return list(cluster_means.var().sort_values(ascending=False).head(n).index)

    sys.exit(f"Unknown --select-regs: {method}")

def aggregate(df, meta, level):
    if level == 'cell':
        return df, meta

    if level == 'animal':
        agg = df.join(meta['animal_id']).groupby('animal_id')[df.columns].mean()
        agg_meta = meta.groupby('animal_id').first().loc[agg.index]
        return agg, agg_meta

    if level == 'cluster':
        grp_col = 'louvain' if 'louvain' in meta.columns else 'ct_louvain'
        agg = df.join(meta[grp_col]).groupby(grp_col)[df.columns].mean()
        agg_meta = pd.DataFrame({grp_col: agg.index}, index=agg.index)
        return agg, agg_meta

    if level == 'animal_region':
        if 'animal_id' not in meta.columns or 'region' not in meta.columns:
            sys.exit("animal_id and region required for animal_region level")
        key = meta['animal_id'].astype(str) + '__' + meta['region'].astype(str)
        tmp = df.copy()
        tmp['_key'] = key.values
        agg = tmp.groupby('_key')[df.columns].mean()
        meta_tmp = meta.copy()
        meta_tmp['_key'] = key.values
        agg_meta = meta_tmp.groupby('_key').first().loc[agg.index]
        return agg, agg_meta

    sys.exit(f"Unknown --level: {level}")

def sort_rows(df, meta, sort_by):
    if not sort_by or sort_by == 'none':
        return df, meta
    if sort_by == 'mean_activity':
        order = df.mean(axis=1).argsort().values
        return df.iloc[order], meta.iloc[order]
    if meta is not None and sort_by in meta.columns:
        col = meta[sort_by]
        if pd.api.types.is_numeric_dtype(col):
            order = np.argsort(col.values.astype(float))
        else:
            order = np.argsort(col.values, kind='stable')
        return df.iloc[order], meta.iloc[order]
    print(f"Warning: sort-by column '{sort_by}' not found, skipping")
    return df, meta

def make_sidebar_colors(values):
    if pd.api.types.is_numeric_dtype(values):
        norm = Normalize(vmin=float(values.min()), vmax=float(values.max()))
        colors = cm.viridis(norm(values.values.astype(float)))
        return colors, 'continuous', cm.viridis, norm, None
    else:
        uniq = sorted(values.unique())
        cmap = cm.get_cmap('tab20', len(uniq))
        cat_map = {v: i for i, v in enumerate(uniq)}
        colors = np.array([cmap(cat_map[v]) for v in values])
        return colors, 'categorical', cmap, None, cat_map

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",        required=True,  help="AUCell scores CSV")
    ap.add_argument("--out",          required=True,  help="Output image path")
    ap.add_argument("--meta",         default="",     help="aucell h5ad for metadata")
    ap.add_argument("--level",        default="cell",
                    choices=["cell", "animal", "cluster", "animal_region"])
    ap.add_argument("--select-regs",  default="variable_between",
                    choices=["variable_between", "variable_within",
                             "variable_cluster", "mean_activity", "variable_cell"])
    ap.add_argument("--regulons",     default="",     help="Comma-separated regulon list")
    ap.add_argument("--n-regs",       type=int, default=30)
    ap.add_argument("--sort-by",      default="",     help="Metadata column or mean_activity")
    ap.add_argument("--color-by",     default="",     help="Metadata column for colored sidebar")
    ap.add_argument("--ytick-cols",   default="",
                    help="Comma-separated metadata columns for y-axis labels e.g. animal_id,age")
    ap.add_argument("--n-cells",      type=int, default=500,
                    help="Cells to sample at cell level (0=all)")
    ap.add_argument("--cluster-regs", action="store_true",
                    help="Cluster regulons by correlation")
    ap.add_argument("--cmap",         default="viridis")
    ap.add_argument("--clip",         type=float, default=0,
                    help="Clip scaled [0-1] values at this level e.g. 0.8 saturates top 20%%")
    ap.add_argument("--seed",         type=int, default=42)
    ap.add_argument("--dpi",          type=int, default=200)
    ap.add_argument("--title",        default="")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"Loading scores: {args.input}")
    df = pd.read_csv(args.input, index_col=0)
    print(f"  shape: {df.shape}")

    meta = None
    if args.meta:
        meta = load_meta(args.meta)
        shared = df.index.intersection(meta.index)
        df, meta = df.loc[shared], meta.loc[shared]
        print(f"  shared cells: {len(shared)}")

    # Select regulons
    regs = select_regulons(df, meta, args.select_regs, args.n_regs, args.regulons)
    reg_method = "specified" if args.regulons.strip() else args.select_regs
    print(f"Selected {len(regs)} regulons via {reg_method}")
    df = df[regs]

    # Aggregate
    df, meta = aggregate(df, meta, args.level)
    print(f"Aggregated shape: {df.shape}")

    # Sample cells
    if args.level == 'cell' and 0 < args.n_cells < len(df):
        idx = rng.choice(df.index.values, size=args.n_cells, replace=False)
        df = df.loc[idx]
        meta = meta.loc[meta.index.isin(idx)].loc[df.index]

    # Sort
    df, meta = sort_rows(df, meta, args.sort_by)

    # Cluster regulons
    if args.cluster_regs and len(regs) > 2:
        Z = linkage(pdist(df.values.T, metric='correlation'), method='average')
        df = df.iloc[:, leaves_list(Z)]
        regs = list(df.columns)
        print("Clustered regulons by correlation")

    # Normalize per regulon to [0,1]
    X = df.values.astype(float)
    col_min = X.min(axis=0, keepdims=True)
    col_max = X.max(axis=0, keepdims=True)
    denom = np.where((col_max - col_min) == 0, 1.0, col_max - col_min)
    Xn = (X - col_min) / denom
    if 0 < args.clip < 1:
        Xn = np.clip(Xn, 0, args.clip) / args.clip

    n_rows, n_cols = Xn.shape

    # Y tick labels
    ytick_cols = [c.strip() for c in args.ytick_cols.split(',') if c.strip()]
    if ytick_cols and meta is not None:
        valid = [c for c in ytick_cols if c in meta.columns]
        def fmt(row_id):
            parts = []
            for c in valid:
                v = meta.loc[row_id, c]
                parts.append(f"{v:.1f}yr" if c == 'age' else str(v))
            return ' | '.join(parts)
        yticklabels = [fmt(i) for i in df.index]
    elif n_rows <= 100:
        yticklabels = list(df.index)
    else:
        yticklabels = None

    # Layout
    has_sidebar = bool(args.color_by and meta is not None
                       and args.color_by in meta.columns)
    ytick_w = 1.8 if yticklabels else 0
    fig_w = max(10, 2.5 + 0.35 * n_cols + ytick_w + (0.4 if has_sidebar else 0))
    fig_h = max(5, min(36, 1.5 + (0.28 * n_rows if args.level != 'cell'
                                   else 0.015 * n_rows)))

    ratios = ([0.03] if has_sidebar else []) + [1]
    fig, axes = plt.subplots(1, len(ratios), figsize=(fig_w, fig_h),
                              gridspec_kw={'width_ratios': ratios, 'wspace': 0.02}
                              if len(ratios) > 1 else None)
    if len(ratios) == 1:
        axes = [axes]

    ax_idx = 0

    # Sidebar
    if has_sidebar:
        ax_side = axes[ax_idx]; ax_idx += 1
        colors, kind, cmap, norm, cat_map = make_sidebar_colors(meta[args.color_by])
        ax_side.imshow(colors.reshape(-1, 1, 4), aspect='auto', interpolation='nearest')
        ax_side.set_xticks([])
        ax_side.set_yticks([])
        ax_side.set_title(args.color_by, fontsize=8, pad=3)
        if kind == 'categorical':
            handles = [Patch(color=cmap(i), label=str(v)) for v, i in cat_map.items()]
            ax_side.legend(handles=handles, bbox_to_anchor=(1.0, 1.0),
                          loc='upper left', fontsize=6, frameon=False)

    # Main heatmap
    ax = axes[ax_idx]
    im = ax.imshow(Xn, aspect='auto', interpolation='nearest', cmap=args.cmap)
    cbar = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.02, shrink=0.5)
    cbar.set_label("AUCell (min-max scaled per regulon;\nabsolute values not comparable)",
                   fontsize=7)

    xfont = max(6, min(10, 200 // n_cols))
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(regs, rotation=55, ha='right', fontsize=xfont)

    if yticklabels:
        yfont = max(5, min(9, 350 // n_rows))
        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(yticklabels, fontsize=yfont)
    else:
        ax.set_yticks([])

    title = args.title or (
        f"AUCell — {n_rows} {args.level}s × {n_cols} regulons"
        + (f" | sorted by {args.sort_by}" if args.sort_by else "")
        + (f" | {reg_method}" if not args.regulons else ""))
    ax.set_title(title, fontsize=10, pad=8)

    plt.savefig(args.out, dpi=args.dpi, bbox_inches='tight')
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
