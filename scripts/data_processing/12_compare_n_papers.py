#!/usr/bin/env python3
"""
Compare n_papers between original and new datasets

Compare the n_papers values between:
- Original: final_balanced_panel_inclusive_by_year_quarter_bi_month_2015_with_groups.dta
- New: final_panel_CLEAN.dta

Focus on gene-months that exist in both datasets to identify differences.
"""

import pandas as pd
import numpy as np

def load_datasets():
    """Load both datasets and standardize column names."""
    print("Loading datasets...")
    
    # Load original dataset
    print("  Loading original dataset...")
    df_orig = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/intermediate_data/final_balanced_panel_inclusive_by_year_quarter_bi_month_2015_with_groups.dta')
    print(f"    Original shape: {df_orig.shape}")
    
    # Load new dataset  
    print("  Loading new dataset...")
    df_new = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/final/final_panel_CLEAN.dta')
    print(f"    New shape: {df_new.shape}")
    
    return df_orig, df_new

def standardize_datasets(df_orig, df_new):
    """Standardize column names and ym formats for comparison."""
    print("Standardizing datasets...")
    
    # Standardize original dataset
    df_orig_std = df_orig.copy()
    df_orig_std = df_orig_std.rename(columns={'Gene_ID': 'gene_id'})
    
    # Convert original ym (datetime) to new format (sequential numbers)
    # Original ym appears to be datetime, new ym is 720-767
    # We need to map them properly
    print("  Converting ym formats...")
    
    # Get unique ym values from original (should be monthly dates)
    orig_ym_unique = sorted(df_orig_std['ym'].unique())
    print(f"    Original ym range: {orig_ym_unique[0]} to {orig_ym_unique[-1]}")
    print(f"    Original unique ym values: {len(orig_ym_unique)}")
    
    # Get unique ym values from new (720-767)
    new_ym_unique = sorted(df_new['ym'].unique())
    print(f"    New ym range: {new_ym_unique[0]:.0f} to {new_ym_unique[-1]:.0f}")
    print(f"    New unique ym values: {len(new_ym_unique)}")
    
    # Create mapping from original datetime to new sequential numbers
    if len(orig_ym_unique) == len(new_ym_unique):
        ym_mapping = dict(zip(orig_ym_unique, new_ym_unique))
        df_orig_std['ym_standardized'] = df_orig_std['ym'].map(ym_mapping)
        print(f"    Successfully mapped {len(ym_mapping)} time periods")
    else:
        print(f"    ⚠️ Warning: Different number of time periods - using month extraction")
        # Fallback: extract year-month and map to sequential
        df_orig_std['year'] = pd.to_datetime(df_orig_std['ym']).dt.year
        df_orig_std['month'] = pd.to_datetime(df_orig_std['ym']).dt.month
        df_orig_std['ym_standardized'] = (df_orig_std['year'] - 2020) * 12 + df_orig_std['month'] + 719
    
    # Use standardized ym
    df_orig_std['ym'] = df_orig_std['ym_standardized']
    
    return df_orig_std, df_new

def find_common_genes(df_orig, df_new):
    """Find genes that exist in both datasets."""
    print("Finding common genes...")
    
    orig_genes = set(df_orig['gene_id'].unique())
    new_genes = set(df_new['gene_id'].unique())
    
    common_genes = orig_genes & new_genes
    orig_only = orig_genes - new_genes
    new_only = new_genes - orig_genes
    
    print(f"  Original genes: {len(orig_genes):,}")
    print(f"  New genes: {len(new_genes):,}")
    print(f"  Common genes: {len(common_genes):,}")
    print(f"  Original only: {len(orig_only):,}")
    print(f"  New only: {len(new_only):,}")
    
    if len(orig_only) > 0:
        print(f"    Sample original-only genes: {list(orig_only)[:5]}")
    if len(new_only) > 0:
        print(f"    Sample new-only genes: {list(new_only)[:5]}")
    
    return common_genes

def compare_n_papers(df_orig, df_new, common_genes):
    """Compare n_papers values for common gene-month combinations."""
    print("Comparing n_papers for common genes...")
    
    # Filter to common genes
    df_orig_common = df_orig[df_orig['gene_id'].isin(common_genes)].copy()
    df_new_common = df_new[df_new['gene_id'].isin(common_genes)].copy()
    
    print(f"  Original records for common genes: {len(df_orig_common):,}")
    print(f"  New records for common genes: {len(df_new_common):,}")
    
    # Merge on gene_id and ym
    merged = df_orig_common.merge(
        df_new_common[['gene_id', 'ym', 'n_papers']], 
        on=['gene_id', 'ym'], 
        how='inner',
        suffixes=('_orig', '_new')
    )
    
    print(f"  Successfully merged records: {len(merged):,}")
    
    if len(merged) == 0:
        print("  ⚠️ No records merged - check ym alignment")
        return None
    
    # Compare n_papers
    merged['n_papers_diff'] = merged['n_papers_new'] - merged['n_papers_orig']
    
    # Statistics
    print(f"\\n=== N_PAPERS COMPARISON ===")
    print(f"Records compared: {len(merged):,}")
    print(f"Identical n_papers: {(merged['n_papers_diff'] == 0).sum():,} ({(merged['n_papers_diff'] == 0).sum() / len(merged) * 100:.1f}%)")
    print(f"Different n_papers: {(merged['n_papers_diff'] != 0).sum():,} ({(merged['n_papers_diff'] != 0).sum() / len(merged) * 100:.1f}%)")
    
    if (merged['n_papers_diff'] != 0).sum() > 0:
        print(f"\\nDifferences summary:")
        print(f"  Mean difference: {merged['n_papers_diff'].mean():.3f}")
        print(f"  Std difference: {merged['n_papers_diff'].std():.3f}")
        print(f"  Min difference: {merged['n_papers_diff'].min()}")
        print(f"  Max difference: {merged['n_papers_diff'].max()}")
        
        # Show distribution of differences
        diff_dist = merged['n_papers_diff'].value_counts().sort_index()
        print(f"\\nTop differences by frequency:")
        print(diff_dist.head(10))
        
        # Show some example differences
        print(f"\\nSample records with differences:")
        diff_samples = merged[merged['n_papers_diff'] != 0].head(10)
        print(diff_samples[['gene_id', 'ym', 'n_papers_orig', 'n_papers_new', 'n_papers_diff']].to_string())
        
        # Check if differences are systematic
        print(f"\\nSystematic patterns:")
        print(f"  Always higher in new: {(merged['n_papers_diff'] > 0).sum():,}")
        print(f"  Always lower in new: {(merged['n_papers_diff'] < 0).sum():,}")
        print(f"  Mixed: {((merged['n_papers_diff'] > 0).sum() > 0) and ((merged['n_papers_diff'] < 0).sum() > 0)}")
    
    return merged

def analyze_missing_records(df_orig, df_new, common_genes):
    """Analyze gene-month combinations that exist in one dataset but not the other."""
    print("\\n=== ANALYZING MISSING RECORDS ===")
    
    # Create keys for comparison
    orig_keys = set(zip(df_orig[df_orig['gene_id'].isin(common_genes)]['gene_id'], 
                       df_orig[df_orig['gene_id'].isin(common_genes)]['ym']))
    new_keys = set(zip(df_new[df_new['gene_id'].isin(common_genes)]['gene_id'], 
                      df_new[df_new['gene_id'].isin(common_genes)]['ym']))
    
    orig_only_keys = orig_keys - new_keys
    new_only_keys = new_keys - orig_keys
    
    print(f"Gene-month combinations in original only: {len(orig_only_keys):,}")
    print(f"Gene-month combinations in new only: {len(new_only_keys):,}")
    
    if len(orig_only_keys) > 0:
        print(f"\\nSample original-only combinations:")
        for i, (gene_id, ym) in enumerate(list(orig_only_keys)[:5]):
            print(f"  Gene {gene_id}, ym {ym}")
    
    if len(new_only_keys) > 0:
        print(f"\\nSample new-only combinations:")
        for i, (gene_id, ym) in enumerate(list(new_only_keys)[:5]):
            print(f"  Gene {gene_id}, ym {ym:.0f}")

def main():
    print("=== N_PAPERS COMPARISON ANALYSIS ===\\n")
    
    # Step 1: Load datasets
    df_orig, df_new = load_datasets()
    
    # Step 2: Standardize formats
    print("\\n" + "="*50)
    df_orig_std, df_new_std = standardize_datasets(df_orig, df_new)
    
    # Step 3: Find common genes
    print("\\n" + "="*50)
    common_genes = find_common_genes(df_orig_std, df_new_std)
    
    # Step 4: Compare n_papers
    print("\\n" + "="*50)
    comparison_result = compare_n_papers(df_orig_std, df_new_std, common_genes)
    
    # Step 5: Analyze missing records
    print("\\n" + "="*50)
    analyze_missing_records(df_orig_std, df_new_std, common_genes)
    
    print(f"\\n=== ANALYSIS COMPLETE ===")

if __name__ == "__main__":
    main()