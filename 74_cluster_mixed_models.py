#!/scratch/easmit31/conda_envs/pyscenic_final/bin/python
"""
Run mixed models separately for each louvain cluster
"""
import sys
import h5py
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def main(in_h5ad, source_h5ad, out_tsv):
    print(f"Loading chromosome info from {source_h5ad}")
    
    with h5py.File(source_h5ad, 'r') as f:
        gene_group = f['var']['external_gene_name']
        gene_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in gene_group['categories'][:]]
        gene_codes = gene_group['codes'][:]
        gene_names = [gene_cats[c] for c in gene_codes]
        
        chr_group = f['var']['chr']
        chr_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in chr_group['categories'][:]]
        chr_codes = chr_group['codes'][:]
        chr_data = [chr_cats[c] for c in chr_codes]
    
    gene_to_chr = dict(zip(gene_names, chr_data))
    
    print(f"\nLoading data from {in_h5ad}")
    
    with h5py.File(in_h5ad, 'r') as f:
        # Get all obs keys that are not standard metadata
        meta_keys = {'age', 'sex', 'animal_id', 'louvain', 'region', 'batch',
                     '_index', 'n_genes', 'n_counts', 'percent_mito'}
        auc_cols = [c for c in f['obs'].keys() if c not in meta_keys]
        
        print(f"\nFound {len(auc_cols)} regulon columns")
        
        # Identify sex chromosome regulons
        print("Identifying sex chromosome regulons...")
        sex_chr_tfs = []
        for reg in auc_cols:
            tf = reg.split('(')[0].split('_extended')[0]
            if tf in gene_to_chr:
                chrom = gene_to_chr[tf]
                if chrom in ['X', 'Y', 'chrX', 'chrY']:
                    sex_chr_tfs.append(tf)
                    print(f"  Excluding {tf} (chr {chrom})")
        
        auc_cols = [c for c in auc_cols if not any(tf in c for tf in sex_chr_tfs)]
        print(f"\nRemoved {len(sex_chr_tfs)} sex chr TFs, keeping {len(auc_cols)} autosomal regulons\n")
        
        # Load metadata
        age = f['obs']['age'][:]
        
        sex_group = f['obs']['sex']
        sex_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in sex_group['categories'][:]]
        sex_codes = sex_group['codes'][:]
        sex = [sex_cats[c] for c in sex_codes]
        
        animal_group = f['obs']['animal_id']
        animal_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in animal_group['categories'][:]]
        animal_codes = animal_group['codes'][:]
        animal_id = [animal_cats[c] for c in animal_codes]
        
        louvain_group = f['obs']['louvain']
        louvain_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in louvain_group['categories'][:]]
        louvain_codes = louvain_group['codes'][:]
        louvain = [louvain_cats[c] for c in louvain_codes]
        
        unique_clusters = sorted(set(louvain))
        print(f"Found {len(unique_clusters)} louvain clusters: {unique_clusters}")
        print(f"Processing {len(auc_cols)} regulons\n")
        
        results = []
        
        for cluster in unique_clusters:
            cluster_mask = [lv == cluster for lv in louvain]
            n_cells = sum(cluster_mask)
            
            print(f"Processing cluster {cluster} ({n_cells} cells)...")
            
            cluster_animals = [animal_id[i] for i, m in enumerate(cluster_mask) if m]
            n_animals = len(set(cluster_animals))
            
            if n_animals < 3:
                print(f"  [skip] Only {n_animals} animals in cluster {cluster}")
                continue
            
            for reg in auc_cols:
                activity = f['obs'][reg][:]
                
                df = pd.DataFrame({
                    "activity": [activity[i] for i, m in enumerate(cluster_mask) if m],
                    "age": [age[i] for i, m in enumerate(cluster_mask) if m],
                    "sex": [sex[i] for i, m in enumerate(cluster_mask) if m],
                    "animal_id": [animal_id[i] for i, m in enumerate(cluster_mask) if m]
                })
                
                df['age'] = pd.to_numeric(df['age'], errors='coerce')
                df['sex'] = pd.Categorical(df['sex'])
                df['animal_id'] = pd.Categorical(df['animal_id'])
                df = df.dropna()
                
                if len(df) < 20:
                    continue
                
                df['age_std'] = (df['age'] - df['age'].mean()) / df['age'].std()
                df['activity_std'] = (df['activity'] - df['activity'].mean()) / df['activity'].std()
                
                try:
                    model = smf.mixedlm(
                        "activity_std ~ age_std + sex",
                        df,
                        groups=df["animal_id"]
                    ).fit(method='powell')
                    
                    re_var = model.cov_re.iloc[0,0] if len(model.cov_re) > 0 else 0
                    resid_var = model.scale
                    icc = re_var / (re_var + resid_var)
                    
                    for term in model.fe_params.index:
                        results.append({
                            "cluster": cluster,
                            "regulon": reg,
                            "term": term,
                            "coef": model.fe_params[term],
                            "se": model.bse[term],
                            "pval": model.pvalues[term],
                            "icc": icc,
                            "n_cells": len(df),
                            "n_animals": df['animal_id'].nunique(),
                            "qval": np.nan
                        })
                        
                except Exception as e:
                    continue
        
        results_df = pd.DataFrame(results)
        
        print("\nApplying FDR correction per term across all clusters...")
        for term in results_df['term'].unique():
            if term == 'Intercept':
                continue
            mask = results_df['term'] == term
            pvals = results_df.loc[mask, 'pval'].values
            if len(pvals) > 0:
                reject, qvals, _, _ = multipletests(pvals, method="fdr_bh")
                results_df.loc[mask, 'qval'] = qvals
                print(f"  {term}: corrected {len(qvals)} p-values")
        
        results_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"\n[done] Saved to {out_tsv}")
        
        print("\nSummary of significant hits (q < 0.05):")
        for term in ['age_std', 'sex[T.M]']:
            sig = results_df[(results_df['term'] == term) & (results_df['qval'] < 0.05)]
            print(f"\n{term}: {len(sig)} significant regulon×cluster combinations")
            if len(sig) > 0:
                print(sig[['cluster', 'regulon', 'coef', 'qval']].head(10))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: ./74_cluster_mixed_models.py aucell.h5ad source.h5ad output.tsv")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
