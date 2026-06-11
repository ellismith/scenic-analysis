#!/scratch/easmit31/conda_envs/pyscenic_final/bin/python
"""
70_regulon_mixed_models.py - FIXED VERSION
"""
import sys
import h5py
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats

def main(in_h5ad, out_tsv):
    print(f"Loading data from {in_h5ad}")
    
    with h5py.File(in_h5ad, 'r') as f:
        meta_keys = {'age','sex','animal_id','louvain','region','batch','_index','n_genes','n_counts','percent_mito'}
        auc_cols = [c for c in f['obs'].keys() if c not in meta_keys and isinstance(f['obs'][c], h5py.Dataset)]
        if not auc_cols:
            raise ValueError("No AUC columns found")
        
        print(f"Found {len(auc_cols)} regulons")
        
        age = f['obs']['age'][:]
        
        sex_group = f['obs']['sex']
        sex_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in sex_group['categories'][:]]
        sex_codes = sex_group['codes'][:]
        sex = [sex_cats[c] for c in sex_codes]
        
        animal_group = f['obs']['animal_id']
        animal_cats = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in animal_group['categories'][:]]
        animal_codes = animal_group['codes'][:]
        animal_id = [animal_cats[c] for c in animal_codes]
        
        results = []
        
        # First, check for age-sex confounding
        df_check = pd.DataFrame({
            "age": pd.to_numeric(age, errors='coerce'),
            "sex": pd.Categorical(sex),
            "animal_id": pd.Categorical(animal_id)
        }).dropna()
        
        print(f"\nData summary:")
        print(f"N observations: {len(df_check)}")
        print(f"N animals: {df_check['animal_id'].nunique()}")
        print(f"\nAge by sex:")
        print(df_check.groupby('sex')['age'].describe())
        
        for reg in auc_cols:
            activity = f['obs'][reg][:]
            
            df = pd.DataFrame({
                "activity": activity,
                "age": pd.to_numeric(age, errors='coerce'),
                "sex": pd.Categorical(sex),
                "animal_id": pd.Categorical(animal_id)
            })
            
            df = df.dropna()
            
            # CRITICAL FIX 1: Standardize age (z-score)
            df['age_std'] = (df['age'] - df['age'].mean()) / df['age'].std()
            
            # CRITICAL FIX 2: Standardize activity within each regulon
            df['activity_std'] = (df['activity'] - df['activity'].mean()) / df['activity'].std()
            
            try:
                # Use standardized variables
                model = smf.mixedlm(
                    "activity_std ~ age_std + sex", 
                    df, 
                    groups=df["animal_id"]
                ).fit(method='powell')  # More robust optimizer
                
                # Extract variance components
                re_var = model.cov_re.iloc[0,0] if len(model.cov_re) > 0 else 0
                resid_var = model.scale
                icc = re_var / (re_var + resid_var)  # Intraclass correlation
                
                for term in model.fe_params.index:
                    results.append({
                        "regulon": reg,
                        "term": term,
                        "coef": model.fe_params[term],
                        "se": model.bse[term],
                        "pval": model.pvalues[term],
                        "random_effect_var": re_var,
                        "residual_var": resid_var,
                        "icc": icc,
                        "converged": model.converged
                    })
                    
            except Exception as e:
                print(f"[warn] {reg} failed: {e}", file=sys.stderr)
                continue
    
    results_df = pd.DataFrame(results)
    
    # Check for convergence issues
    if 'converged' in results_df.columns:
        n_unconverged = (~results_df['converged']).sum()
        if n_unconverged > 0:
            print(f"WARNING: {n_unconverged} models did not converge!")
    
    # Check ICC distribution
    if 'icc' in results_df.columns:
        print(f"\nICC summary (proportion of variance from random effects):")
        print(results_df.groupby('term')['icc'].describe())
    
    results_df['qval'] = np.nan
    
    # Apply FDR correction separately for age and sex, but NOT intercept
    print("\nApplying FDR correction per term (excluding intercept)...")
    for term in results_df['term'].unique():
        # Skip intercept from FDR correction
        if term == 'Intercept':
            continue
            
        mask = results_df['term'] == term
        pvals = results_df.loc[mask, 'pval'].values
        reject, qvals, alphacSidak, alphacBonf = multipletests(pvals, method="fdr_bh")
        results_df.loc[mask, 'qval'] = qvals
    
    results_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[done] Saved mixed-effects results to {out_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: ./70_regulon_mixed_models.py input.h5ad output.tsv")
    main(sys.argv[1], sys.argv[2])
