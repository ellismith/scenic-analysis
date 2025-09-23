#!/usr/bin/env python
import os, argparse, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

p = argparse.ArgumentParser(description="TF modules via Jaccard overlap of targets from adj.csv")
p.add_argument("--adj", default="adj.csv")
p.add_argument("--min_weight", type=float, default=0.0, help="Keep edges with importance >= this (default 0)")
p.add_argument("--top_tf", type=int, default=60, help="Use top-N TFs by target count (default 60)")
p.add_argument("--dist_thresh", type=float, default=0.7, help="Distance threshold (1-Jaccard) to cut clusters")
p.add_argument("--out_fig", default="figs/adj_tf_modules_jaccard.png")
p.add_argument("--out_tab_dir", default="tables")
args = p.parse_args()
os.makedirs(os.path.dirname(args.out_fig), exist_ok=True)
os.makedirs(args.out_tab_dir, exist_ok=True)

df = pd.read_csv(args.adj)
cols = {c.lower(): c for c in df.columns}
TFc = cols.get("tf") or cols.get("regulator") or list(df.columns)[0]
TGc = cols.get("target") or cols.get("gene") or list(df.columns)[1]
IMc = cols.get("importance") or cols.get("weight") or list(df.columns)[2]

# filter & build TF->targets (unique)
if IMc in df:
    df = df[df[IMc].astype(float) >= args.min_weight]
tf_to_targets = (
    df.groupby(TFc)[TGc]
      .apply(lambda s: set(map(str, s.values)))
      .to_dict()
)

# pick top-N TFs by target count
counts = {tf: len(tg) for tf, tg in tf_to_targets.items()}
top = sorted(counts.keys(), key=lambda t: counts[t], reverse=True)[:args.top_tf]
top = [t for t in top if counts[t] > 0]
if len(top) < 2:
    raise SystemExit("Not enough TFs with targets after filtering.")

# Jaccard similarity matrix
n = len(top)
sets = [tf_to_targets[t] for t in top]
J = np.zeros((n, n), dtype=float)
for i in range(n):
    Ai = sets[i]
    for j in range(i, n):
        Aj = sets[j]
        inter = len(Ai & Aj)
        union = len(Ai | Aj)
        J[i, j] = J[j, i] = (inter / union) if union else 0.0

# cluster on distance = 1 - J
D = 1.0 - J
Z = linkage(squareform(D, checks=False), method="average")
order = leaves_list(Z)
Jr = J[order][:, order]
ordered_tfs = [top[i] for i in order]

# cut into modules
labels = fcluster(Z, t=args.dist_thresh, criterion="distance")
labels = labels[order]  # reorder labels to match heatmap
out_clusters = pd.DataFrame({"TF": ordered_tfs, "cluster": labels, "n_targets": [counts[t] for t in ordered_tfs]})
out_clusters.to_csv(os.path.join(args.out_tab_dir, "tf_modules_clusters.tsv"), sep="\t", index=False)

# save TF target counts
pd.DataFrame({"TF": list(counts.keys()), "n_targets": list(counts.values())}) \
  .sort_values("n_targets", ascending=False) \
  .to_csv(os.path.join(args.out_tab_dir, "tf_target_counts.tsv"), sep="\t", index=False)

# plot heatmap
fig_w = max(6, 0.16 * n + 3)
fig_h = max(6, 0.16 * n + 2)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
im = ax.imshow(Jr, vmin=0, vmax=1, aspect="equal")
ax.set_title(f"TF modules (Jaccard of target sets)\nN TFs={n}, min_weight≥{args.min_weight}")
# tick every ~5 to avoid clutter
step = max(1, n // 30)
ax.set_xticks(range(0, n, step)); ax.set_xticklabels([ordered_tfs[i] for i in range(0, n, step)], rotation=90, fontsize=7)
ax.set_yticks(range(0, n, step)); ax.set_yticklabels([ordered_tfs[i] for i in range(0, n, step)], fontsize=7)
cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("Jaccard similarity")
plt.tight_layout()
fig.savefig(args.out_fig, dpi=200)
plt.close(fig)

print("[OK] wrote:", args.out_fig)
print("[OK] clusters ->", os.path.join(args.out_tab_dir, "tf_modules_clusters.tsv"))
print("[OK] counts   ->", os.path.join(args.out_tab_dir, "tf_target_counts.tsv"))
