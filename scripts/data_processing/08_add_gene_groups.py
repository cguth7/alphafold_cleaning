#!/usr/bin/env python3
"""
Phase 8: Add Gene Groups to Final Panel

Merge gene group classifications onto the balanced panel.
Gene groups provide functional/pathway classifications for genes.

Input: final_panel_COMPLETE.parquet + genes_groups.dta
Output: final_panel_WITH_GROUPS.parquet/.dta
"""

import pandas as pd
import numpy as np
from pathlib import Path

def main():
    print("=== PHASE 8: ADDING GENE GROUPS ===\n")
    
    # Load final balanced panel
    print("Loading final balanced panel...")
    df_panel = pd.read_parquet('final/final_panel_COMPLETE.parquet')
    print(f"  Panel: {len(df_panel):,} observations")
    print(f"  Unique genes: {df_panel['gene_id'].nunique():,}")
    print(f"  Gene ID type: {df_panel['gene_id'].dtype}")
    print(f"  Gene ID range: {df_panel['gene_id'].min()}-{df_panel['gene_id'].max()}")
    
    # Load gene groups
    print("\nLoading gene groups data...")
    df_groups = pd.read_stata('intermediate_data/genes_groups.dta')
    print(f"  Gene groups: {len(df_groups):,} records")
    print(f"  Columns: {df_groups.columns.tolist()}")
    print(f"  Gene ID type: {df_groups['ncbi_entrez_gene_id'].dtype}")
    print(f"  Gene ID range: {df_groups['ncbi_entrez_gene_id'].min()}-{df_groups['ncbi_entrez_gene_id'].max()}")
    print(f"  Unique genes with groups: {df_groups['ncbi_entrez_gene_id'].nunique():,}")
    
    # Check for missing/empty gene groups
    print(f"\nGene group data quality:")
    missing_groups = df_groups['GenegroupID'].isna().sum()
    empty_groups = (df_groups['GenegroupID'] == '').sum()
    print(f"  Missing gene groups: {missing_groups:,}")
    print(f"  Empty gene groups: {empty_groups:,}")
    print(f"  Valid gene groups: {len(df_groups) - missing_groups - empty_groups:,}")
    
    # Sample of gene groups
    print(f"\nSample gene groups:")
    valid_groups = df_groups[df_groups['GenegroupID'].notna() & (df_groups['GenegroupID'] != '')]
    if len(valid_groups) > 0:
        print(valid_groups['GenegroupID'].value_counts().head(10))
    
    # Convert gene_id types for merging
    print(f"\nPreparing merge...")
    df_groups['gene_id'] = df_groups['ncbi_entrez_gene_id'].astype('int64')
    
    # Check overlap between panel and gene groups
    panel_genes = set(df_panel['gene_id'].unique())
    group_genes = set(df_groups['gene_id'].unique())
    overlap = panel_genes & group_genes
    
    print(f"  Panel genes: {len(panel_genes):,}")
    print(f"  Gene group genes: {len(group_genes):,}")
    print(f"  Overlap: {len(overlap):,} genes ({len(overlap)/len(panel_genes)*100:.1f}% of panel)")
    
    # Merge gene groups onto panel
    print(f"\nMerging gene groups...")
    merged = df_panel.merge(
        df_groups[['gene_id', 'GenegroupID']],
        on='gene_id',
        how='left'
    )
    
    print(f"  Merged panel: {len(merged):,} observations")
    
    # Check merge success
    genes_with_groups = merged['GenegroupID'].notna().sum()
    genes_with_valid_groups = merged[
        merged['GenegroupID'].notna() & (merged['GenegroupID'] != '')
    ].shape[0]
    
    print(f"  Observations with gene groups: {genes_with_groups:,} ({genes_with_groups/len(merged)*100:.1f}%)")
    print(f"  Observations with valid gene groups: {genes_with_valid_groups:,} ({genes_with_valid_groups/len(merged)*100:.1f}%)")
    
    # Gene-level coverage
    unique_genes_with_groups = merged.drop_duplicates('gene_id')['GenegroupID'].notna().sum()
    unique_genes_with_valid_groups = merged.drop_duplicates('gene_id')[
        merged.drop_duplicates('gene_id')['GenegroupID'].notna() & 
        (merged.drop_duplicates('gene_id')['GenegroupID'] != '')
    ].shape[0]
    
    print(f"  Unique genes with groups: {unique_genes_with_groups:,} ({unique_genes_with_groups/df_panel['gene_id'].nunique()*100:.1f}%)")
    print(f"  Unique genes with valid groups: {unique_genes_with_valid_groups:,} ({unique_genes_with_valid_groups/df_panel['gene_id'].nunique()*100:.1f}%)")
    
    # Handle multiple gene groups (some entries have formats like "442|454")
    print(f"\nHandling multiple gene groups...")
    multi_group_pattern = merged['GenegroupID'].str.contains('\\|', na=False)
    multi_groups = multi_group_pattern.sum()
    print(f"  Observations with multiple groups (|): {multi_groups:,}")
    
    if multi_groups > 0:
        print("  Sample multiple groups:")
        sample_multi = merged[multi_group_pattern]['GenegroupID'].value_counts().head(5)
        print(sample_multi)
    
    # Rename column to be more descriptive
    merged = merged.rename(columns={'GenegroupID': 'gene_group'})
    
    # Final validation
    print(f"\n=== FINAL VALIDATION ===")
    print(f"Panel structure maintained: {len(merged) == len(df_panel)} ({'✅' if len(merged) == len(df_panel) else '❌'})")
    print(f"All original columns preserved: {set(df_panel.columns).issubset(set(merged.columns))} ({'✅' if set(df_panel.columns).issubset(set(merged.columns)) else '❌'})")
    print(f"Gene group coverage: {unique_genes_with_valid_groups:,}/{df_panel['gene_id'].nunique():,} genes ({unique_genes_with_valid_groups/df_panel['gene_id'].nunique()*100:.1f}%)")
    
    # Show final column list
    print(f"\nFinal dataset columns ({len(merged.columns)}):")
    for i, col in enumerate(merged.columns):
        print(f"  {i+1:2d}. {col}")
    
    # Save enhanced dataset
    print(f"\n=== SAVING ENHANCED DATASET ===")
    
    # Parquet version
    parquet_path = Path("final/final_panel_WITH_GROUPS.parquet")
    merged.to_parquet(parquet_path, index=False)
    parquet_size = parquet_path.stat().st_size / 1024**2
    print(f"Saved: {parquet_path} ({parquet_size:.1f} MB)")
    
    # Also update the COMPLETE version with gene groups
    complete_path = Path("final/final_panel_COMPLETE.parquet")
    merged.to_parquet(complete_path, index=False)
    print(f"Updated: {complete_path} (now includes gene_group column)")
    
    # Summary for documentation
    print(f"\n🎉 SUCCESS! Gene groups added to final panel!")
    print(f"   - Total observations: {len(merged):,} (balanced)")
    print(f"   - Gene group coverage: {unique_genes_with_valid_groups:,}/{df_panel['gene_id'].nunique():,} genes ({unique_genes_with_valid_groups/df_panel['gene_id'].nunique()*100:.1f}%)")
    print(f"   - Multiple group assignments: {multi_groups:,} observations")
    print(f"   - Ready for functional analysis of AlphaFold impact!")

if __name__ == "__main__":
    main()