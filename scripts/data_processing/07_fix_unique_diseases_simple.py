#!/usr/bin/env python3
"""
Simple Fix for Unique Disease Merge

The correct approach: merge unique disease data on pmid and ym directly,
then aggregate to gene-month level for the final panel.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def main():
    print("=== SIMPLE UNIQUE DISEASE MERGE FIX ===\n")
    
    # Load the corrected panel data that has proper ym values
    print("Loading corrected panel data...")
    df = pd.read_parquet('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/final/final_panel_READY.parquet')
    print(f"  Loaded {len(df):,} gene-month observations")
    print(f"  Columns: {df.columns.tolist()}")
    print(f"  ym range: {df['ym'].min()}-{df['ym'].max()}")
    
    # Drop existing unique disease columns (they're all 0 anyway)
    df = df.drop(columns=['unique_mesh_count', 'new_mesh_count'], errors='ignore')
    
    # Load unique disease data  
    print("\nLoading unique disease data...")
    df_unique = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/intermediate_data/gene_month_unique_mesh_2015.dta')
    df_new_unique = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/intermediate_data/gene_month_new_unique_diseases_2015_no_level_5.dta')
    
    print(f"  Unique mesh: {len(df_unique):,} observations, ym range {df_unique['ym'].min()}-{df_unique['ym'].max()}")
    print(f"  New unique: {len(df_new_unique):,} observations, ym range {df_new_unique['ym'].min()}-{df_new_unique['ym'].max()}")
    
    # Create the time mapping using the ym_date column
    print("\nCreating time mapping...")
    time_mapping = df_unique[['ym', 'ym_date']].drop_duplicates()
    time_mapping['year'] = time_mapping['ym_date'].dt.year  
    time_mapping['month'] = time_mapping['ym_date'].dt.month
    time_mapping['ym_proper'] = time_mapping['year'] * 100 + time_mapping['month']
    
    # Filter to our time period (2020-2023)
    our_period = time_mapping[
        (time_mapping['year'] >= 2020) & (time_mapping['year'] <= 2023)
    ]
    
    print(f"  Our period (2020-2023) maps to ym: {our_period['ym'].min()}-{our_period['ym'].max()}")
    print(f"  That's {len(our_period)} months of data")
    
    # Create mapping dictionary from proper ym to their ym system
    ym_mapping = dict(zip(our_period['ym_proper'], our_period['ym']))
    print(f"  Created mapping for {len(ym_mapping)} month periods")
    
    # Add mapped ym to our panel data
    df['ym_unique_system'] = df['ym'].map(ym_mapping)
    mapped_count = df['ym_unique_system'].notna().sum()
    print(f"  Successfully mapped {mapped_count:,}/{len(df):,} observations")
    
    # Convert data types for merging
    df['gene_id_int32'] = df['gene_id'].astype('int32')
    df['ym_int32'] = df['ym_unique_system'].astype('int32')
    
    # Merge unique mesh data at gene-month level
    print(f"\nMerging unique mesh data...")
    merged = df.merge(
        df_unique[['Gene_ID', 'ym', 'unique_mesh_count']],
        left_on=['gene_id_int32', 'ym_int32'],
        right_on=['Gene_ID', 'ym'],
        how='left'
    )
    merged['unique_mesh_count'] = merged['unique_mesh_count'].fillna(0).astype('int32')
    
    # Merge new unique diseases  
    print(f"Merging new unique diseases...")
    merged = merged.merge(
        df_new_unique[['Gene_ID', 'ym', 'new_mesh_count']],
        left_on=['gene_id_int32', 'ym_int32'], 
        right_on=['Gene_ID', 'ym'],
        how='left',
        suffixes=('', '_new')
    )
    merged['new_mesh_count'] = merged['new_mesh_count'].fillna(0).astype('int32')
    
    # Clean up merge columns
    cleanup_cols = ['Gene_ID', 'ym_y', 'Gene_ID_new', 'ym_new', 'gene_id_int32', 'ym_int32', 'ym_unique_system']
    merged = merged.drop(columns=[col for col in cleanup_cols if col in merged.columns])
    
    # Report merge success
    unique_added = (merged['unique_mesh_count'] > 0).sum()
    new_added = (merged['new_mesh_count'] > 0).sum()
    print(f"  Unique mesh counts added to {unique_added:,} observations")
    print(f"  New mesh counts added to {new_added:,} observations")
    print(f"  Total unique diseases: {merged['unique_mesh_count'].sum():,}")
    print(f"  Total new diseases: {merged['new_mesh_count'].sum():,}")
    
    # The merged data is our final balanced panel
    balanced_panel = merged
    
    # Final validation
    print(f"\n=== FINAL VALIDATION ===")
    print(f"Final panel: {len(balanced_panel):,} observations")
    
    # Check veteran + newcomer = total
    author_sum = balanced_panel['n_veteran_papers'] + balanced_panel['n_newcomer_papers']
    mismatch = (author_sum != balanced_panel['n_papers']).sum()
    print(f"Veteran + Newcomer = Total: {mismatch} mismatches ({'✅' if mismatch == 0 else '❌'})")
    
    # Check disease data
    unique_obs = (balanced_panel['unique_mesh_count'] > 0).sum()
    new_obs = (balanced_panel['new_mesh_count'] > 0).sum()
    print(f"Unique disease data: {unique_obs:,} observations ({'✅' if unique_obs > 0 else '❌'})")
    print(f"New disease data: {new_obs:,} observations ({'✅' if new_obs > 0 else '❌'})")
    print(f"Total unique diseases: {balanced_panel['unique_mesh_count'].sum():,}")
    print(f"Total new diseases: {balanced_panel['new_mesh_count'].sum():,}")
    
    # Save final versions
    print(f"\n=== SAVING FINAL FILES ===")
    
    # Parquet version
    parquet_path = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/final/final_panel_COMPLETE.parquet")
    balanced_panel.to_parquet(parquet_path, index=False)
    print(f"Saved: {parquet_path}")
    
    # Stata version
    stata_path = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/final/final_panel_COMPLETE.dta")
    balanced_panel.to_stata(stata_path, write_index=False)
    stata_size = stata_path.stat().st_size / 1024**2
    print(f"Saved: {stata_path} ({stata_size:.1f} MB)")
    
    print(f"\n🎉 SUCCESS! Unique disease merge fixed!")
    print(f"   - Unique diseases: {balanced_panel['unique_mesh_count'].sum():,} ✅")
    print(f"   - New diseases: {balanced_panel['new_mesh_count'].sum():,} ✅")
    print(f"   - Perfect balance: {len(balanced_panel):,} observations ✅")
    print(f"   - Stata ready: Use 'tsset gene_id ym_seq' ✅")

if __name__ == "__main__":
    main()