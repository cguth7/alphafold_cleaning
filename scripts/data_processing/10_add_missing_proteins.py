#!/usr/bin/env python3
"""
Phase 10: Add Missing Proteins from Master Data

Add the 361 missing proteins from protein_data_master to the final panel.
These proteins will be added with a dummy flag to distinguish them from the main dataset.

Input:  
- final/final_panel_WITH_DEPOSITS.dta (945K rows, current panel)
- intermediate_data/protein_data_master_dta.dta (1.2M rows, master data)

Output: final/final_panel_COMPLETE_WITH_MISSING.dta and .parquet

Processing steps:
1. Load current final panel and identify existing proteins
2. Load master data and identify missing proteins  
3. Create panel structure for missing proteins (2020-2023)
4. Add dummy flag and fill appropriate values
5. Combine with existing panel
6. Save complete dataset
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_current_panel(file_path):
    """Load the current final panel."""
    print(f"Loading current final panel from: {file_path}")
    
    df = pd.read_stata(file_path)
    
    print(f"  Loaded {len(df):,} panel records")
    print(f"  Unique proteins: {df['protein_id'].nunique():,}")
    print(f"  Time span: {df['year'].min()}-{df['year'].max()}")
    print(f"  Columns: {len(df.columns)}")
    
    return df

def identify_missing_proteins(df_panel, df_master):
    """Identify proteins in master that are missing from panel."""
    print("Identifying missing proteins...")
    
    panel_proteins = set(df_panel['protein_id'].unique())
    master_proteins = set(df_master['protein_id'].unique())
    missing_proteins = master_proteins - panel_proteins
    
    print(f"  Panel proteins: {len(panel_proteins):,}")
    print(f"  Master proteins: {len(master_proteins):,}")
    print(f"  Missing proteins: {len(missing_proteins):,}")
    
    if len(missing_proteins) > 0:
        print(f"  Sample missing: {list(missing_proteins)[:5]}")
    
    return missing_proteins

def create_missing_protein_panel(missing_proteins, df_master, df_panel_template):
    """Create panel structure for missing proteins."""
    print(f"Creating panel structure for {len(missing_proteins):,} missing proteins...")
    
    # Get the time structure from existing panel (2020-2023, 48 months)
    # Filter to only complete records (some have NaN in ym column)
    time_structure = df_panel_template[['year', 'month', 'quarter', 'year_quarter', 'bimonth', 'year_bimonth', 'ym_x', 'ym_seq', 'ym']].dropna().drop_duplicates().sort_values('ym_x')
    print(f"  Time periods: {len(time_structure)} months")
    
    # Get protein info for missing proteins
    missing_protein_info = df_master[df_master['protein_id'].isin(missing_proteins)].copy()
    
    # Get unique protein characteristics (take first occurrence for each protein)
    protein_chars = missing_protein_info.groupby('protein_id').agg({
        'average_plddt': 'first',
        'protein_existence': 'first', 
        'gene_name': 'first',
        'geneid': 'first'
    }).reset_index()
    
    print(f"  Missing protein characteristics collected: {len(protein_chars):,}")
    
    # Create full panel structure (protein × time)
    missing_panel_records = []
    
    for _, time_row in time_structure.iterrows():
        for _, protein_row in protein_chars.iterrows():
            record = {
                # Time variables
                'year': time_row['year'],
                'month': time_row['month'],
                'quarter': time_row['quarter'],
                'year_quarter': time_row['year_quarter'],
                'bimonth': time_row['bimonth'],
                'year_bimonth': time_row['year_bimonth'],
                'ym_x': time_row['ym_x'],
                'ym_seq': time_row['ym_seq'],
                'ym': time_row['ym'],
                
                # Protein characteristics
                'protein_id': protein_row['protein_id'],
                'gene_id': protein_row['geneid'],
                'average_plddt': protein_row['average_plddt'],
                'protein_existence': protein_row['protein_existence'],
                'gene_name': protein_row['gene_name'],
                
                # Publication metrics (all zeros - no publications found)
                'n_papers': 0,
                'n_top25_y': 0,
                'n_top10_y': 0,
                'n_top05_y': 0,
                'n_top25_q': 0,
                'n_top10_q': 0,
                'n_top05_q': 0,
                'n_top25_b2': 0,
                'n_top10_b2': 0,
                'n_top05_b2': 0,
                'n_newcomer_papers': 0,
                'n_veteran_papers': 0,
                
                # Disease metrics (all zeros)
                'unique_mesh_count': 0,
                'new_mesh_count': 0,
                
                # Deposits data (all zeros - assume no deposits)
                'num_deposits': 0,
                'monthly_deposits': 0,
                
                # Gene group (missing)
                'gene_group': None,
                
                # Dummy flag
                'missing_from_literature': 1
            }
            missing_panel_records.append(record)
    
    df_missing = pd.DataFrame(missing_panel_records)
    
    expected_records = len(missing_proteins) * len(time_structure)
    print(f"  Created {len(df_missing):,} records ({expected_records:,} expected)")
    print(f"  Missing proteins × {len(time_structure)} months = {len(missing_proteins)} × {len(time_structure)} = {expected_records}")
    
    return df_missing

def add_dummy_flag_to_existing(df_panel):
    """Add dummy flag to existing panel data."""
    print("Adding dummy flag to existing panel...")
    
    df_with_flag = df_panel.copy()
    df_with_flag['missing_from_literature'] = 0  # Not missing
    
    print(f"  Added missing_from_literature flag to {len(df_with_flag):,} existing records")
    
    return df_with_flag

def combine_panels(df_existing, df_missing):
    """Combine existing panel with missing protein panel."""
    print("Combining existing and missing protein panels...")
    
    # Ensure columns match
    existing_cols = set(df_existing.columns)
    missing_cols = set(df_missing.columns)
    
    print(f"  Existing panel columns: {len(existing_cols)}")
    print(f"  Missing panel columns: {len(missing_cols)}")
    
    if existing_cols != missing_cols:
        print(f"  Column mismatch detected:")
        only_existing = existing_cols - missing_cols
        only_missing = missing_cols - existing_cols
        if only_existing:
            print(f"    Only in existing: {only_existing}")
        if only_missing:
            print(f"    Only in missing: {only_missing}")
        
        # Align columns
        all_cols = existing_cols | missing_cols
        for col in all_cols:
            if col not in df_existing.columns:
                df_existing[col] = None
            if col not in df_missing.columns:
                df_missing[col] = None
    
    # Combine
    df_combined = pd.concat([df_existing, df_missing], ignore_index=True)
    
    print(f"  Combined panel: {len(df_combined):,} records")
    print(f"  Existing records: {(df_combined['missing_from_literature'] == 0).sum():,}")
    print(f"  Missing protein records: {(df_combined['missing_from_literature'] == 1).sum():,}")
    print(f"  Total proteins: {df_combined['protein_id'].nunique():,}")
    
    return df_combined

def validate_complete_panel(df):
    """Validate the complete panel structure."""
    print("Validating complete panel...")
    
    total_records = len(df)
    total_proteins = df['protein_id'].nunique()
    time_periods = df['ym_x'].nunique()
    
    expected_records = total_proteins * time_periods
    
    print(f"  Total records: {total_records:,}")
    print(f"  Total proteins: {total_proteins:,}")
    print(f"  Time periods: {time_periods}")
    print(f"  Expected records: {expected_records:,}")
    print(f"  Panel balanced: {total_records == expected_records}")
    
    # Check missing flag distribution
    missing_records = (df['missing_from_literature'] == 1).sum()
    existing_records = (df['missing_from_literature'] == 0).sum()
    
    print(f"  Missing protein records: {missing_records:,} ({missing_records/total_records*100:.1f}%)")
    print(f"  Existing records: {existing_records:,} ({existing_records/total_records*100:.1f}%)")
    
    # Validate that missing proteins have zero publications
    missing_subset = df[df['missing_from_literature'] == 1]
    has_papers = (missing_subset['n_papers'] > 0).sum()
    print(f"  Missing proteins with papers: {has_papers} (should be 0)")
    
    return total_records == expected_records

def save_complete_panel(df, base_path):
    """Save the complete panel with missing proteins."""
    print(f"Saving complete panel...")
    
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
    summary_path = base_path.parent / 'missing_proteins_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("Missing Proteins Addition Summary\\n")
        f.write("=" * 35 + "\\n\\n")
        
        f.write(f"Process: Added {df['protein_id'].nunique() - (df['missing_from_literature'] == 0).sum() // 48} missing proteins from master data\\n\\n")
        
        f.write(f"Input files:\\n")
        f.write(f"  - final/final_panel_WITH_DEPOSITS.dta (existing panel)\\n")
        f.write(f"  - intermediate_data/protein_data_master_dta.dta (master data)\\n\\n")
        
        f.write(f"Output files:\\n")
        f.write(f"  - {dta_path.name} ({dta_size:.1f} MB)\\n")
        f.write(f"  - {parquet_path.name} ({parquet_size:.1f} MB)\\n\\n")
        
        f.write(f"Final dataset summary:\\n")
        f.write(f"  Total records: {len(df):,}\\n")
        f.write(f"  Total proteins: {df['protein_id'].nunique():,}\\n")
        f.write(f"  Time periods: 48 months (2020-2023)\\n")
        f.write(f"  Records from literature: {(df['missing_from_literature'] == 0).sum():,}\\n")
        f.write(f"  Records for missing proteins: {(df['missing_from_literature'] == 1).sum():,}\\n\\n")
        
        f.write(f"New flag added:\\n")
        f.write(f"  - missing_from_literature: 1 for proteins not found in literature, 0 for existing\\n")
    
    print(f"  Summary saved: {summary_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    current_panel_input = base_dir / "final" / "final_panel_WITH_DEPOSITS.dta"
    master_data_input = base_dir / "intermediate_data" / "protein_data_master_dta.dta"
    output_base = base_dir / "final" / "final_panel_COMPLETE_WITH_MISSING"
    
    print("=== Phase 10: Add Missing Proteins ===\\n")
    
    # Step 1: Load current panel
    df_panel = load_current_panel(current_panel_input)
    
    # Step 2: Load master data
    print("\\n" + "="*50)
    print("Loading master protein data...")
    df_master = pd.read_stata(master_data_input)
    print(f"  Loaded {len(df_master):,} master records")
    print(f"  Unique proteins in master: {df_master['protein_id'].nunique():,}")
    
    # Step 3: Identify missing proteins
    print("\\n" + "="*50)
    missing_proteins = identify_missing_proteins(df_panel, df_master)
    
    if len(missing_proteins) == 0:
        print("No missing proteins found. Current panel is complete.")
        return
    
    # Step 4: Create panel structure for missing proteins
    print("\\n" + "="*50)
    df_missing = create_missing_protein_panel(missing_proteins, df_master, df_panel)
    
    # Step 5: Add dummy flag to existing panel
    print("\\n" + "="*50)
    df_existing_with_flag = add_dummy_flag_to_existing(df_panel)
    
    # Step 6: Combine panels
    print("\\n" + "="*50)
    df_complete = combine_panels(df_existing_with_flag, df_missing)
    
    # Step 7: Validate complete panel
    print("\\n" + "="*50)
    is_valid = validate_complete_panel(df_complete)
    
    if not is_valid:
        print("⚠️ Warning: Panel validation failed. Check data integrity.")
    
    # Step 8: Save complete panel
    print("\\n" + "="*50)
    save_complete_panel(df_complete, output_base)
    
    print(f"\\n=== Phase 10 Complete ===")
    print(f"Complete panel: {len(df_complete):,} records")
    print(f"Total proteins: {df_complete['protein_id'].nunique():,}")
    print(f"Missing proteins added: {len(missing_proteins):,}")
    print(f"Output: final_panel_COMPLETE_WITH_MISSING.dta and .parquet")

if __name__ == "__main__":
    main()