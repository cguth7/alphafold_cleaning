#!/usr/bin/env python3
"""
Phase 9: Merge Deposits Data into Final Panel

Merge deposits data from PDB into the final panel:
1. deposits_dta.dta - merge on protein_id (total deposits per protein)
2. deposits_monthly_dta.dta - merge on protein_id and ym (monthly deposits)

Input:  
- final/final_panel_WITH_GROUPS.dta (945K rows, final panel)
- intermediate_data/deposits_dta.dta (7.8K rows, total deposits)
- intermediate_data/deposits_monthly_dta.dta (46.7K rows, monthly deposits)

Output: final/final_panel_WITH_DEPOSITS.dta and .parquet

Processing steps:
1. Load final panel data
2. Load and examine deposits data structures
3. Merge deposits_dta on protein_id 
4. Convert monthly deposits dates to ym format and merge
5. Handle missing values and validate merges
6. Save enhanced final panel
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_final_panel(file_path):
    """Load the current final panel."""
    print(f"Loading final panel from: {file_path}")
    
    df = pd.read_stata(file_path)
    
    print(f"  Loaded {len(df):,} panel records")
    print(f"  Panel dimensions: {df.shape}")
    print(f"  Unique proteins: {df['protein_id'].nunique():,}")
    print(f"  Time span: {df['year'].min()}-{df['year'].max()}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def load_deposits_data(total_deposits_path, monthly_deposits_path):
    """Load both deposits datasets."""
    print(f"Loading deposits data...")
    
    # Load total deposits
    print(f"  Loading total deposits from: {total_deposits_path}")
    df_total = pd.read_stata(total_deposits_path)
    print(f"    Shape: {df_total.shape}")
    print(f"    Columns: {list(df_total.columns)}")
    print(f"    Sample data:")
    print(f"    {df_total.head(3).to_string()}")
    
    # Load monthly deposits
    print(f"  Loading monthly deposits from: {monthly_deposits_path}")
    df_monthly = pd.read_stata(monthly_deposits_path)
    print(f"    Shape: {df_monthly.shape}")
    print(f"    Columns: {list(df_monthly.columns)}")
    print(f"    Date range: {df_monthly['month_num'].min()} to {df_monthly['month_num'].max()}")
    print(f"    Sample data:")
    print(f"    {df_monthly.head(3).to_string()}")
    
    return df_total, df_monthly

def prepare_monthly_deposits(df_monthly):
    """Convert monthly deposits to ym format for merging."""
    print("Preparing monthly deposits for merging...")
    
    # Convert month_num to ym format (YYYYMM)
    df_monthly['ym_deposits'] = df_monthly['month_num'].dt.year * 100 + df_monthly['month_num'].dt.month
    
    print(f"  Converted dates to ym format:")
    print(f"    ym range: {df_monthly['ym_deposits'].min()} - {df_monthly['ym_deposits'].max()}")
    
    # Filter to our panel time period (2020-2023)
    panel_period = df_monthly['ym_deposits'].between(202001, 202312)
    df_filtered = df_monthly[panel_period].copy()
    
    print(f"  Filtered to panel period (2020-2023):")
    print(f"    Records: {len(df_filtered):,} ({len(df_filtered)/len(df_monthly)*100:.1f}% of monthly data)")
    print(f"    Proteins: {df_filtered['protein_id'].nunique():,}")
    
    # Prepare for merge (rename columns to avoid conflicts)
    df_filtered = df_filtered[['protein_id', 'ym_deposits', 'num_deposits_month']].rename(columns={
        'ym_deposits': 'ym_x',  # Match the format in final panel
        'num_deposits_month': 'monthly_deposits'
    })
    
    return df_filtered

def merge_total_deposits(df_panel, df_deposits):
    """Merge total deposits data with panel."""
    print("Merging total deposits data...")
    
    initial_records = len(df_panel)
    initial_proteins = df_panel['protein_id'].nunique()
    
    # Check overlap
    panel_proteins = set(df_panel['protein_id'].unique())
    deposits_proteins = set(df_deposits['protein_id'].unique())
    overlap = panel_proteins & deposits_proteins
    
    print(f"  Panel proteins: {len(panel_proteins):,}")
    print(f"  Deposits proteins: {len(deposits_proteins):,}")
    print(f"  Overlap: {len(overlap):,} proteins ({len(overlap)/len(panel_proteins)*100:.1f}% of panel)")
    
    # Merge
    df_merged = df_panel.merge(df_deposits[['protein_id', 'num_deposits']], on='protein_id', how='left')
    
    # Fill missing values with 0 (no deposits)
    df_merged['num_deposits'] = df_merged['num_deposits'].fillna(0)
    
    final_records = len(df_merged)
    deposits_coverage = (df_merged['num_deposits'] > 0).sum()
    
    print(f"  Merged records: {final_records:,} (should match panel: {final_records == initial_records})")
    print(f"  Records with deposits: {deposits_coverage:,} ({deposits_coverage/final_records*100:.1f}%)")
    print(f"  Average deposits per protein: {df_merged['num_deposits'].mean():.2f}")
    print(f"  Max deposits: {df_merged['num_deposits'].max():.0f}")
    
    return df_merged

def merge_monthly_deposits(df_panel, df_monthly):
    """Merge monthly deposits data with panel."""
    print("Merging monthly deposits data...")
    
    initial_records = len(df_panel)
    
    # Check overlap in both protein_id and ym_x
    panel_keys = set(zip(df_panel['protein_id'], df_panel['ym_x']))
    monthly_keys = set(zip(df_monthly['protein_id'], df_monthly['ym_x']))
    overlap = panel_keys & monthly_keys
    
    print(f"  Panel protein-month combinations: {len(panel_keys):,}")
    print(f"  Monthly deposits combinations: {len(monthly_keys):,}")
    print(f"  Overlap: {len(overlap):,} combinations ({len(overlap)/len(panel_keys)*100:.1f}% of panel)")
    
    # Merge
    df_merged = df_panel.merge(df_monthly, on=['protein_id', 'ym_x'], how='left')
    
    # Fill missing values with 0 (no monthly deposits)
    df_merged['monthly_deposits'] = df_merged['monthly_deposits'].fillna(0)
    
    final_records = len(df_merged)
    monthly_coverage = (df_merged['monthly_deposits'] > 0).sum()
    
    print(f"  Merged records: {final_records:,} (should match panel: {final_records == initial_records})")
    print(f"  Records with monthly deposits: {monthly_coverage:,} ({monthly_coverage/final_records*100:.1f}%)")
    print(f"  Average monthly deposits: {df_merged['monthly_deposits'].mean():.3f}")
    print(f"  Max monthly deposits: {df_merged['monthly_deposits'].max():.0f}")
    
    return df_merged

def validate_merged_data(df):
    """Validate the merged data quality."""
    print("Validating merged data...")
    
    # Check data integrity
    total_records = len(df)
    print(f"  Total records: {total_records:,}")
    print(f"  Unique proteins: {df['protein_id'].nunique():,}")
    print(f"  Time span: {df['year'].min()}-{df['year'].max()}")
    
    # Check deposits data
    has_total_deposits = (df['num_deposits'] > 0).sum()
    has_monthly_deposits = (df['monthly_deposits'] > 0).sum()
    
    print(f"  Records with total deposits: {has_total_deposits:,} ({has_total_deposits/total_records*100:.1f}%)")
    print(f"  Records with monthly deposits: {has_monthly_deposits:,} ({has_monthly_deposits/total_records*100:.1f}%)")
    
    # Summary statistics
    print(f"  Deposits summary:")
    print(f"    Total deposits: mean={df['num_deposits'].mean():.2f}, median={df['num_deposits'].median():.0f}")
    print(f"    Monthly deposits: mean={df['monthly_deposits'].mean():.3f}, median={df['monthly_deposits'].median():.0f}")
    
    # Check for any correlation
    correlation = df['num_deposits'].corr(df.groupby('protein_id')['monthly_deposits'].sum())
    print(f"  Correlation between total and summed monthly deposits: {correlation:.3f}")
    
    return True

def save_enhanced_panel(df, base_path):
    """Save the enhanced panel with deposits data."""
    print(f"Saving enhanced panel...")
    
    # Save as both .dta and .parquet
    dta_path = base_path.with_suffix('.dta')
    parquet_path = base_path.with_suffix('.parquet')
    
    print(f"  Saving Stata format: {dta_path}")
    df.to_stata(dta_path, write_index=False)
    dta_size = os.path.getsize(dta_path) / 1024**2
    
    print(f"  Saving Parquet format: {parquet_path}")
    df.to_parquet(parquet_path, index=False)
    parquet_size = os.path.getsize(parquet_path) / 1024**2
    
    print(f"  File sizes:")
    print(f"    {dta_path.name}: {dta_size:.1f} MB")
    print(f"    {parquet_path.name}: {parquet_size:.1f} MB")
    
    # Save processing summary
    summary_path = base_path.parent / 'deposits_merge_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("Deposits Data Merge Summary\\n")
        f.write("=" * 30 + "\\n\\n")
        
        f.write(f"Input files:\\n")
        f.write(f"  - final/final_panel_WITH_GROUPS.dta (original panel)\\n")
        f.write(f"  - intermediate_data/deposits_dta.dta (total deposits)\\n")
        f.write(f"  - intermediate_data/deposits_monthly_dta.dta (monthly deposits)\\n\\n")
        
        f.write(f"Output files:\\n")
        f.write(f"  - {dta_path.name} ({dta_size:.1f} MB)\\n")
        f.write(f"  - {parquet_path.name} ({parquet_size:.1f} MB)\\n\\n")
        
        f.write(f"Data summary:\\n")
        f.write(f"  Total records: {len(df):,}\\n")
        f.write(f"  Proteins with total deposits: {(df['num_deposits'] > 0).sum():,}\\n")
        f.write(f"  Records with monthly deposits: {(df['monthly_deposits'] > 0).sum():,}\\n\\n")
        
        f.write(f"New columns added:\\n")
        f.write(f"  - num_deposits: Total deposits per protein\\n")
        f.write(f"  - monthly_deposits: Deposits per protein per month\\n")
    
    print(f"  Summary saved: {summary_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    final_panel_input = base_dir / "final" / "final_panel_WITH_GROUPS.dta"
    total_deposits_input = base_dir / "intermediate_data" / "deposits_dta.dta"
    monthly_deposits_input = base_dir / "intermediate_data" / "deposits_monthly_dta.dta"
    output_base = base_dir / "final" / "final_panel_WITH_DEPOSITS"
    
    print("=== Phase 9: Merge Deposits Data ===\\n")
    
    # Step 1: Load final panel
    df_panel = load_final_panel(final_panel_input)
    
    # Step 2: Load deposits data
    print("\\n" + "="*50)
    df_total_deposits, df_monthly_deposits = load_deposits_data(total_deposits_input, monthly_deposits_input)
    
    # Step 3: Prepare monthly deposits for merging
    print("\\n" + "="*50)
    df_monthly_prepared = prepare_monthly_deposits(df_monthly_deposits)
    
    # Step 4: Merge total deposits
    print("\\n" + "="*50)
    df_with_total = merge_total_deposits(df_panel, df_total_deposits)
    
    # Step 5: Merge monthly deposits
    print("\\n" + "="*50)
    df_final = merge_monthly_deposits(df_with_total, df_monthly_prepared)
    
    # Step 6: Validate results
    print("\\n" + "="*50)
    validate_merged_data(df_final)
    
    # Step 7: Save enhanced panel
    print("\\n" + "="*50)
    save_enhanced_panel(df_final, output_base)
    
    print(f"\\n=== Phase 9 Complete ===")
    print(f"Enhanced panel: {len(df_final):,} records")
    print(f"New deposits columns: num_deposits, monthly_deposits")
    print(f"Output: final_panel_WITH_DEPOSITS.dta and .parquet")

if __name__ == "__main__":
    main()