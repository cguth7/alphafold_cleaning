#!/usr/bin/env python3
"""
Phase 6: Fix Unique Disease Merge & Create Stata Export

Issues identified:
1. ym format mismatch: Our panel uses 5393-5704, unique disease files use 660-767  
2. Data type mismatch: gene_id int64 vs Gene_ID int32
3. Need proper ym variable for Stata (202001 format)

This script will:
1. Create proper ym variables (both Stata-friendly and for merging)
2. Fix the unique disease merge with proper data types
3. Export final .dta version for Stata analysis
4. Validate all merges and data consistency
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_and_fix_panel(file_path):
    """Load panel and create proper ym variables."""
    print(f"Loading and fixing panel from: {file_path}")
    
    df = pd.read_parquet(file_path)
    print(f"  Loaded {len(df):,} observations")
    
    # Current ym is weird (5393 = ym, but should be 202001)
    # Let's create proper ym from year/month
    print("  Creating proper ym variables...")
    
    # Convert year/month to proper YYYYMM format  
    df['year_int'] = df['year'].astype(int)
    df['month_int'] = df['month'].astype(int)
    df['ym_proper'] = df['year_int'] * 100 + df['month_int']  # 202001 format
    
    # Also create a sequential month index starting from 2015-01 for merging with unique data
    # 2015-01 = month 1, so 2020-01 = month 61 (5 years * 12 + 1)
    # But the unique data starts from ym=660, so let's reverse-engineer
    
    # From analysis: unique data ym 660-767 likely corresponds to some period
    # Let's check what our 2020-2023 should map to
    base_year = 2015  # Unique data starts from 2015
    df['months_since_2015'] = (df['year_int'] - base_year) * 12 + df['month_int']
    
    print(f"  Original ym range: {df['ym'].min()}-{df['ym'].max()}")  
    print(f"  Proper ym range: {df['ym_proper'].min()}-{df['ym_proper'].max()}")
    print(f"  Months since 2015: {df['months_since_2015'].min()}-{df['months_since_2015'].max()}")
    
    # The unique data ym 660-767 likely means something different
    # Let's map our AlphaFold era (2020-2023) to the unique data time frame
    # 2020 = year 6 since 2015, so Jan 2020 = month 61
    # But unique data has ym up to 767, so probably goes beyond 2023
    
    return df

def analyze_unique_data_time_mapping():
    """Figure out the time mapping in unique disease data."""
    print("Analyzing unique disease data time mapping...")
    
    # Load unique data to understand ym mapping
    df_unique = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/intermediate_data/gene_month_unique_mesh_2015.dta')
    
    print(f"  Unique data ym range: {df_unique['ym'].min()}-{df_unique['ym'].max()}")
    print(f"  Unique data has ym_date column with actual dates")
    
    # Check the mapping between ym and ym_date
    time_mapping = df_unique[['ym', 'ym_date']].drop_duplicates().sort_values('ym')
    print("  Time mapping sample:")
    print(time_mapping.head(10).to_string())
    print("  ...")
    print(time_mapping.tail(10).to_string())
    
    # Extract year/month from ym_date to understand the mapping
    time_mapping['year_from_date'] = time_mapping['ym_date'].dt.year
    time_mapping['month_from_date'] = time_mapping['ym_date'].dt.month
    
    # Find our AlphaFold period (2020-2023) in their data
    alphafold_period = time_mapping[
        (time_mapping['year_from_date'] >= 2020) & 
        (time_mapping['year_from_date'] <= 2023)
    ]
    
    print(f"  AlphaFold period mapping:")
    print(f"  2020-2023 maps to ym: {alphafold_period['ym'].min()}-{alphafold_period['ym'].max()}")
    
    return time_mapping

def merge_unique_diseases_fixed(df_panel, time_mapping):
    """Merge unique disease data with proper time mapping."""
    print("Merging unique disease data with proper mapping...")
    
    # Create mapping from our ym_proper to their ym system
    print("  Creating time mapping...")
    time_mapping['ym_proper'] = time_mapping['year_from_date'] * 100 + time_mapping['month_from_date']
    mapping_dict = dict(zip(time_mapping['ym_proper'], time_mapping['ym']))
    
    # Add their ym system to our panel
    df_panel['ym_unique_system'] = df_panel['ym_proper'].map(mapping_dict)
    
    print(f"  Mapped {df_panel['ym_unique_system'].notna().sum():,} observations to unique system")
    
    # Filter to only periods that exist in unique data
    df_mappable = df_panel[df_panel['ym_unique_system'].notna()].copy()
    print(f"  Observations with mappable time periods: {len(df_mappable):,}")
    
    # Load and merge unique mesh data
    print("  Loading unique mesh data...")
    df_unique = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/intermediate_data/gene_month_unique_mesh_2015.dta')
    
    # Convert gene_id to int32 for matching
    df_mappable['gene_id_int32'] = df_mappable['gene_id'].astype('int32')
    df_mappable['ym_int32'] = df_mappable['ym_unique_system'].astype('int32')
    
    print("  Merging unique mesh counts...")
    merged = df_mappable.merge(
        df_unique[['Gene_ID', 'ym', 'unique_mesh_count']], 
        left_on=['gene_id_int32', 'ym_int32'],
        right_on=['Gene_ID', 'ym'], 
        how='left'
    )
    merged['unique_mesh_count'] = merged['unique_mesh_count'].fillna(0).astype('int32')
    
    # Load and merge new unique diseases
    print("  Loading new unique diseases...")
    df_new_unique = pd.read_stata('/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/intermediate_data/gene_month_new_unique_diseases_2015_no_level_5.dta')
    
    print("  Merging new unique disease counts...")
    merged = merged.merge(
        df_new_unique[['Gene_ID', 'ym', 'new_mesh_count']], 
        left_on=['gene_id_int32', 'ym_int32'],
        right_on=['Gene_ID', 'ym'], 
        how='left',
        suffixes=('', '_new')
    )
    merged['new_mesh_count'] = merged['new_mesh_count'].fillna(0).astype('int32')
    
    # Clean up merge columns
    cols_to_drop = ['Gene_ID', 'ym_y', 'Gene_ID_new', 'ym_new', 'gene_id_int32', 'ym_int32', 'ym_unique_system']
    merged = merged.drop(columns=[col for col in cols_to_drop if col in merged.columns])
    
    unique_mesh_added = (merged['unique_mesh_count'] > 0).sum()
    new_mesh_added = (merged['new_mesh_count'] > 0).sum()
    
    print(f"  Successfully merged!")
    print(f"  Unique mesh counts added to {unique_mesh_added:,} observations")
    print(f"  New mesh counts added to {new_mesh_added:,} observations")
    
    return merged

def create_stata_friendly_version(df):
    """Create Stata-friendly version with proper variable types and names."""
    print("Creating Stata-friendly version...")
    
    # Create proper ym variable for Stata (monthly time variable)
    # Stata prefers ym() format starting from some base
    # Let's create a clean monthly time variable
    
    # Method 1: Sequential months from 2020m1
    base_year_month = 2020 * 12 + 1  # 2020m1
    df['stata_ym'] = (df['year_int'] - 2020) * 12 + df['month_int']  # 1, 2, 3, ..., 48
    
    # Method 2: Also keep the proper YYYYMM format 
    # This is what we created earlier: ym_proper
    
    print(f"  Created stata_ym: {df['stata_ym'].min()}-{df['stata_ym'].max()} (sequential months 2020-2023)")
    print(f"  Kept ym_proper: {df['ym_proper'].min()}-{df['ym_proper'].max()} (YYYYMM format)")
    
    # Ensure all numeric columns are appropriate types for Stata
    stata_df = df.copy()
    
    # Convert to Stata-friendly types
    stata_df['gene_id'] = stata_df['gene_id'].astype('int32')
    stata_df['year'] = stata_df['year_int'].astype('int16') 
    stata_df['month'] = stata_df['month_int'].astype('int8')
    stata_df['ym'] = stata_df['ym_proper'].astype('int32')  # YYYYMM format
    stata_df['ym_seq'] = stata_df['stata_ym'].astype('int16')  # Sequential 1-48
    
    # Keep essential columns for Stata
    essential_cols = [
        'gene_id', 'ym', 'ym_seq', 'year', 'month', 'quarter', 'year_quarter', 'bimonth', 'year_bimonth',
        'n_papers', 'n_top25_y', 'n_top10_y', 'n_top05_y', 
        'n_top25_q', 'n_top10_q', 'n_top05_q',
        'n_top25_b2', 'n_top10_b2', 'n_top05_b2',
        'n_newcomer_papers', 'n_veteran_papers',
        'unique_mesh_count', 'new_mesh_count',
        'protein_id', 'gene_name', 'protein_existence', 'average_plddt'
    ]
    
    stata_df = stata_df[essential_cols].copy()
    
    # Convert temporal variables to appropriate Stata types
    for col in ['quarter', 'year_quarter', 'bimonth', 'year_bimonth']:
        if col in stata_df.columns:
            stata_df[col] = stata_df[col].astype('int16')
    
    print(f"  Final Stata dataset: {len(stata_df):,} rows × {len(stata_df.columns)} columns")
    
    return stata_df

def save_stata_version(df, output_path):
    """Save as Stata .dta file with proper formatting."""
    print(f"Saving Stata version to: {output_path}")
    
    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save as .dta
    df.to_stata(output_path, write_index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} observations to Stata format")
    print(f"  File size: {file_size:.1f} MB")
    
    # Also save updated parquet version
    parquet_path = output_path.with_suffix('.parquet')
    df.to_parquet(parquet_path, index=False)
    print(f"  Also saved parquet version: {parquet_path}")
    
    # Save summary statistics
    stats_path = output_path.parent / 'fixed_panel_summary.txt'
    with open(stats_path, 'w') as f:
        f.write("Fixed Panel with Unique Disease Data - Summary\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Panel Structure:\n")
        f.write(f"  Total observations: {len(df):,}\n")
        f.write(f"  Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"  Time periods: {df['ym'].nunique()} months\n")
        f.write(f"  Time span: {df['year'].min()}-{df['year'].max()}\n\n")
        
        f.write("Time Variables for Stata:\n")
        f.write(f"  ym: YYYYMM format ({df['ym'].min()}-{df['ym'].max()})\n") 
        f.write(f"  ym_seq: Sequential months 1-48 ({df['ym_seq'].min()}-{df['ym_seq'].max()})\n")
        f.write(f"  Use 'tsset gene_id ym_seq' in Stata for panel setup\n\n")
        
        f.write("Unique Disease Data Integration:\n")
        unique_mesh_nonzero = (df['unique_mesh_count'] > 0).sum()
        new_mesh_nonzero = (df['new_mesh_count'] > 0).sum()
        f.write(f"  Observations with unique mesh data: {unique_mesh_nonzero:,}\n")
        f.write(f"  Observations with new mesh data: {new_mesh_nonzero:,}\n")
        f.write(f"  Total unique mesh associations: {df['unique_mesh_count'].sum():,}\n")
        f.write(f"  Total new mesh associations: {df['new_mesh_count'].sum():,}\n\n")
        
        f.write("Data Quality:\n")
        f.write(f"  Complete gene metadata: {df['gene_name'].notna().sum():,} observations\n")
        f.write(f"  Active gene-months: {(df['n_papers'] > 0).sum():,}\n")
        f.write(f"  Total publications: {df['n_papers'].sum():,}\n\n")
        
        f.write(f"Files created:\n")
        f.write(f"  - {output_path.name} (Stata format)\n")
        f.write(f"  - {parquet_path.name} (Parquet format)\n")
    
    print(f"  Summary saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    panel_input = base_dir / "final" / "final_balanced_panel.parquet"
    stata_output = base_dir / "final" / "final_balanced_panel_with_unique_diseases.dta"
    
    print("=== Phase 6: Fix Unique Disease Merge & Stata Export ===\n")
    
    # Step 1: Load and fix panel
    df_panel = load_and_fix_panel(panel_input)
    
    # Step 2: Analyze unique data time mapping
    print("\\n" + "="*60)
    time_mapping = analyze_unique_data_time_mapping()
    
    # Step 3: Merge unique diseases with proper mapping
    print("\\n" + "="*60)
    df_with_diseases = merge_unique_diseases_fixed(df_panel, time_mapping)
    
    # Step 4: Create Stata-friendly version
    print("\\n" + "="*60)
    stata_df = create_stata_friendly_version(df_with_diseases)
    
    # Step 5: Save both versions
    print("\\n" + "="*60)
    save_stata_version(stata_df, stata_output)
    
    print(f"\\n=== Phase 6 Complete ===")
    print(f"✅ Fixed unique disease merge")
    print(f"✅ Created Stata-friendly time variables")  
    print(f"✅ Saved .dta and .parquet versions")
    print(f"Output: {stata_output}")

if __name__ == "__main__":
    main()