#!/usr/bin/env python3
"""
Phase 11: Create Clean Final Panel

Create a cleaned version of the final panel that:
1. Drops redundant time columns, keeping only ym
2. Ensures ym is non-null for all records
3. Optimizes column order and data types
4. Saves as .dta format for analysis

Input:  
- final/final_panel_COMPLETE_WITH_MISSING.dta (complete panel)

Output: final/final_panel_CLEAN.dta (cleaned for analysis)

Processing steps:
1. Load complete panel
2. Check ym column completeness
3. Drop redundant time columns
4. Optimize column order and data types
5. Validate cleaned panel
6. Save cleaned version
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_complete_panel(file_path):
    """Load the complete panel with all data."""
    print(f"Loading complete panel from: {file_path}")
    
    df = pd.read_stata(file_path)
    
    print(f"  Loaded {len(df):,} panel records")
    print(f"  Panel dimensions: {df.shape}")
    print(f"  Unique proteins: {df['protein_id'].nunique():,}")
    print(f"  Columns: {list(df.columns)}")
    
    return df

def check_ym_completeness(df):
    """Check that ym column is complete and valid."""
    print("Checking ym column completeness...")
    
    total_records = len(df)
    null_ym = df['ym'].isnull().sum()
    
    print(f"  Total records: {total_records:,}")
    print(f"  Null ym values: {null_ym:,}")
    print(f"  ym completeness: {(total_records - null_ym) / total_records * 100:.2f}%")
    
    if null_ym > 0:
        print(f"  ⚠️ Warning: {null_ym:,} records have null ym values")
        
        # Show sample of null records
        null_records = df[df['ym'].isnull()][['protein_id', 'year', 'month', 'ym_x', 'ym']]
        print(f"  Sample null records:")
        print(null_records.head().to_string())
        
        # Check if we can fill from ym_x
        can_fill = df[df['ym'].isnull() & df['ym_x'].notnull()].shape[0]
        print(f"  Records with null ym but valid ym_x: {can_fill:,}")
        
        return False
    else:
        print(f"  ✓ All records have valid ym values")
        
        # Show ym range
        print(f"  ym range: {df['ym'].min():.0f} - {df['ym'].max():.0f}")
        unique_ym = df['ym'].nunique()
        print(f"  Unique ym values: {unique_ym}")
        
        return True

def fix_missing_ym(df):
    """Fix missing ym values using ym_x if possible."""
    print("Fixing missing ym values...")
    
    null_ym = df['ym'].isnull().sum()
    if null_ym == 0:
        print("  No missing ym values to fix")
        return df
    
    # Try to fill from ym_x
    df_fixed = df.copy()
    mask = df_fixed['ym'].isnull() & df_fixed['ym_x'].notnull()
    
    if mask.sum() > 0:
        # Convert ym_x (YYYYMM) to the ym format used in the dataset
        # ym is a sequential counter: 720 = 202001, 721 = 202002, etc.
        
        # Find the mapping from existing data
        mapping_data = df_fixed[df_fixed['ym'].notnull() & df_fixed['ym_x'].notnull()]
        if len(mapping_data) > 0:
            # Create a mapping dictionary
            ym_mapping = dict(zip(mapping_data['ym_x'], mapping_data['ym']))
            print(f"  Created mapping from {len(ym_mapping)} ym_x values")
            
            # Apply mapping to missing values
            df_fixed.loc[mask, 'ym'] = df_fixed.loc[mask, 'ym_x'].map(ym_mapping)
        else:
            # Fallback: assume 720 = 202001 and increment
            base_ym_x = 202001
            base_ym = 720
            df_fixed.loc[mask, 'ym'] = df_fixed.loc[mask, 'ym_x'] - base_ym_x + base_ym
        
        filled = mask.sum()
        print(f"  Filled {filled:,} missing ym values from ym_x")
        
        # Check if any still missing
        still_missing = df_fixed['ym'].isnull().sum()
        if still_missing > 0:
            print(f"  ⚠️ Still {still_missing:,} records with missing ym")
            # Remove these records
            df_fixed = df_fixed[df_fixed['ym'].notnull()].copy()
            print(f"  Removed records with missing ym. Final count: {len(df_fixed):,}")
    
    return df_fixed

def select_and_order_columns(df):
    """Select essential columns and order them logically."""
    print("Selecting and ordering columns...")
    
    # Define column groups
    core_columns = [
        'gene_id',
        'protein_id', 
        'gene_name',
        'ym'  # Only time column we keep
    ]
    
    protein_info = [
        'average_plddt',
        'protein_existence',
        'gene_group'
    ]
    
    publication_metrics = [
        'n_papers',
        'n_top25_y', 'n_top10_y', 'n_top05_y',
        'n_top25_q', 'n_top10_q', 'n_top05_q', 
        'n_top25_b2', 'n_top10_b2', 'n_top05_b2',
        'n_newcomer_papers',
        'n_veteran_papers'
    ]
    
    disease_metrics = [
        'unique_mesh_count',
        'new_mesh_count'
    ]
    
    deposits_metrics = [
        'num_deposits',
        'monthly_deposits'
    ]
    
    flags = [
        'missing_from_literature'
    ]
    
    # Combine in logical order
    selected_columns = core_columns + protein_info + publication_metrics + disease_metrics + deposits_metrics + flags
    
    # Check which columns exist in the dataframe
    available_columns = []
    missing_columns = []
    
    for col in selected_columns:
        if col in df.columns:
            available_columns.append(col)
        else:
            missing_columns.append(col)
    
    print(f"  Selected {len(available_columns)} available columns")
    if missing_columns:
        print(f"  Missing columns (will skip): {missing_columns}")
    
    # Select available columns
    df_selected = df[available_columns].copy()
    
    print(f"  Final dimensions: {df_selected.shape}")
    print(f"  Column order: {list(df_selected.columns)}")
    
    return df_selected

def optimize_data_types(df):
    """Optimize data types for efficient storage."""
    print("Optimizing data types...")
    
    initial_memory = df.memory_usage(deep=True).sum() / 1024**2
    
    df_optimized = df.copy()
    
    # Optimize integer columns
    int_columns = df_optimized.select_dtypes(include=['int32', 'int64']).columns
    for col in int_columns:
        if col not in ['gene_id']:  # Keep gene_id as-is for safety
            # Check if we can downcast
            col_min = df_optimized[col].min()
            col_max = df_optimized[col].max()
            
            if col_min >= 0:  # Unsigned integers
                if col_max <= 255:
                    df_optimized[col] = df_optimized[col].astype('uint8')
                elif col_max <= 65535:
                    df_optimized[col] = df_optimized[col].astype('uint16')
                elif col_max <= 4294967295:
                    df_optimized[col] = df_optimized[col].astype('uint32')
            else:  # Signed integers
                if col_min >= -128 and col_max <= 127:
                    df_optimized[col] = df_optimized[col].astype('int8')
                elif col_min >= -32768 and col_max <= 32767:
                    df_optimized[col] = df_optimized[col].astype('int16')
    
    # Optimize float columns
    float_columns = df_optimized.select_dtypes(include=['float64']).columns
    for col in float_columns:
        if col != 'ym':  # Keep ym precision for safety
            # Try float32 if precision allows
            df_optimized[col] = pd.to_numeric(df_optimized[col], downcast='float')
    
    final_memory = df_optimized.memory_usage(deep=True).sum() / 1024**2
    memory_reduction = (initial_memory - final_memory) / initial_memory * 100
    
    print(f"  Memory usage: {initial_memory:.1f} MB → {final_memory:.1f} MB ({memory_reduction:.1f}% reduction)")
    
    return df_optimized

def validate_cleaned_panel(df):
    """Validate the cleaned panel structure."""
    print("Validating cleaned panel...")
    
    total_records = len(df)
    total_proteins = df['protein_id'].nunique()
    total_ym = df['ym'].nunique()
    
    print(f"  Total records: {total_records:,}")
    print(f"  Total proteins: {total_proteins:,}")
    print(f"  Total time periods: {total_ym}")
    print(f"  Expected records: {total_proteins * total_ym:,}")
    
    # Check panel balance
    protein_counts = df.groupby('protein_id').size()
    expected_periods = protein_counts.mode()[0] if len(protein_counts) > 0 else 0
    balanced_proteins = (protein_counts == expected_periods).sum()
    
    print(f"  Expected periods per protein: {expected_periods}")
    print(f"  Proteins with expected periods: {balanced_proteins:,} ({balanced_proteins/total_proteins*100:.1f}%)")
    
    # Check ym completeness
    null_ym = df['ym'].isnull().sum()
    print(f"  Records with null ym: {null_ym:,}")
    
    # Check missing_from_literature flag
    if 'missing_from_literature' in df.columns:
        missing_count = (df['missing_from_literature'] == 1).sum()
        existing_count = (df['missing_from_literature'] == 0).sum()
        print(f"  Missing protein records: {missing_count:,}")
        print(f"  Existing protein records: {existing_count:,}")
    
    # Data quality checks
    print(f"  Data quality checks:")
    print(f"    Unique genes: {df['gene_id'].nunique():,}")
    print(f"    Records with publications: {(df['n_papers'] > 0).sum():,}")
    if 'num_deposits' in df.columns:
        print(f"    Records with deposits: {(df['num_deposits'] > 0).sum():,}")
    
    is_valid = (null_ym == 0) and (total_records == total_proteins * total_ym)
    return is_valid

def save_cleaned_panel(df, output_path):
    """Save the cleaned panel."""
    print(f"Saving cleaned panel to: {output_path}")
    
    # Save as .dta format for analysis
    df.to_stata(output_path, write_index=False)
    file_size = os.path.getsize(output_path) / 1024**2
    
    print(f"  File size: {file_size:.1f} MB")
    print(f"  Records: {len(df):,}")
    print(f"  Columns: {len(df.columns)}")
    
    # Save column documentation
    doc_path = output_path.parent / 'clean_panel_columns.txt'
    with open(doc_path, 'w') as f:
        f.write("Clean Final Panel - Column Documentation\\n")
        f.write("=" * 45 + "\\n\\n")
        
        f.write(f"File: {output_path.name}\\n")
        f.write(f"Records: {len(df):,}\\n")
        f.write(f"Proteins: {df['protein_id'].nunique():,}\\n")
        f.write(f"Time periods: {df['ym'].nunique()}\\n\\n")
        
        f.write("Columns retained:\\n")
        for i, col in enumerate(df.columns, 1):
            dtype = str(df[col].dtype)
            non_null = df[col].notna().sum()
            f.write(f"  {i:2d}. {col:<25} ({dtype:<10}) - {non_null:,} non-null\\n")
        
        f.write("\\nColumns removed:\\n")
        removed_cols = ['year', 'month', 'quarter', 'year_quarter', 'bimonth', 'year_bimonth', 'ym_x', 'ym_seq']
        for col in removed_cols:
            f.write(f"  - {col} (redundant time variable)\\n")
        
        f.write("\\nTime variable retained:\\n")
        f.write(f"  - ym: Primary time identifier ({df['ym'].min():.0f} to {df['ym'].max():.0f})\\n")
    
    print(f"  Documentation saved: {doc_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    input_panel = base_dir / "final" / "final_panel_COMPLETE_WITH_MISSING.dta"
    output_panel = base_dir / "final" / "final_panel_CLEAN.dta"
    
    print("=== Phase 11: Create Clean Final Panel ===\\n")
    
    # Step 1: Load complete panel
    df = load_complete_panel(input_panel)
    
    # Step 2: Check ym completeness
    print("\\n" + "="*50)
    ym_complete = check_ym_completeness(df)
    
    # Step 3: Fix missing ym if needed
    if not ym_complete:
        print("\\n" + "="*50)
        df = fix_missing_ym(df)
    
    # Step 4: Select and order columns
    print("\\n" + "="*50)
    df = select_and_order_columns(df)
    
    # Step 5: Optimize data types
    print("\\n" + "="*50)
    df = optimize_data_types(df)
    
    # Step 6: Validate cleaned panel
    print("\\n" + "="*50)
    is_valid = validate_cleaned_panel(df)
    
    if not is_valid:
        print("⚠️ Warning: Cleaned panel validation failed")
    
    # Step 7: Save cleaned panel
    print("\\n" + "="*50)
    save_cleaned_panel(df, output_panel)
    
    print(f"\\n=== Phase 11 Complete ===")
    print(f"Clean panel: {len(df):,} records")
    print(f"Columns: {len(df.columns)} (dropped redundant time variables)")
    print(f"Time variable: ym only")
    print(f"Output: {output_panel.name}")

if __name__ == "__main__":
    main()