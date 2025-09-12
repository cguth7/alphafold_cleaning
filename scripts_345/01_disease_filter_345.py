#!/usr/bin/env python3
"""
Phase 1: Disease-Relevant Publication Filter (3_4_5 Levels)

Filter disease PubTator data to keep only PMIDs with:
1. MESH diseases at depth levels 3/4/5 OR  
2. OMIM diseases

This matches the Stata pipeline exactly, including level 5.

Input:  cleaned/disease_pmid_mesh.txt (163M rows, 3.4GB)  
Output: processed_345/disease_relevant_pmids_345.parquet

Processing steps:
1. Load cleaned disease data (pmid, mesh_id)
2. Strip MESH:/OMIM: prefixes and flag OMIM entries
3. Filter by MESH depth 3/4/5 reference list (NOT just 3/4)
4. Create unique PMID list
5. Save with full documentation
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_disease_data(file_path):
    """Load cleaned disease PubTator data."""
    print(f"Loading disease data from: {file_path}")
    
    # Read tab-separated file with headers
    df = pd.read_csv(file_path, sep='\t', dtype={'pmid': 'int32', 'mesh_id': 'str'})
    
    print(f"  Loaded {len(df):,} rows")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def process_mesh_identifiers(df):
    """Strip prefixes and flag OMIM entries."""
    print("Processing MESH/OMIM identifiers...")
    
    # Flag OMIM entries before stripping prefixes
    df['is_omim'] = df['mesh_id'].str.contains('OMIM:', na=False)
    omim_count = df['is_omim'].sum()
    print(f"  Found {omim_count:,} OMIM entries ({omim_count/len(df)*100:.2f}%)")
    
    # Strip prefixes
    df['mesh_id_clean'] = (df['mesh_id']
                          .str.replace('OMIM:', '', regex=False)
                          .str.replace('MESH:', '', regex=False))
    
    return df

def load_mesh_depth_filter(filter_file):
    """Load the MESH depth 3/4/5 filter list."""
    print(f"Loading MESH depth 3/4/5 filter from: {filter_file}")
    
    # Try to load from the Dropbox path mentioned in Stata
    mesh_file = "/Users/maxguthmann/Dropbox/Matteo_Danilo_AlphaFold2/rawdata/charlie/mesh_depth_3_4_5_no_explode.dta"
    
    if os.path.exists(mesh_file):
        print(f"  Using Stata file: {mesh_file}")
        mesh_df = pd.read_stata(mesh_file)
        
        # Check column names
        print(f"  Columns in mesh file: {list(mesh_df.columns)}")
        
        # Likely column names: DescriptorUI, Mesh_ID, mesh_id, or similar
        mesh_col = None
        for col in ['DescriptorUI', 'Mesh_ID', 'mesh_id', 'mesh_ui', 'descriptor_ui']:
            if col in mesh_df.columns:
                mesh_col = col
                break
        
        if mesh_col is None:
            # Take first column if no obvious match
            mesh_col = mesh_df.columns[0]
            print(f"  Warning: Using first column '{mesh_col}' as mesh identifier")
        else:
            print(f"  Using column '{mesh_col}' as mesh identifier")
        
        # Extract unique mesh IDs
        valid_mesh_ids = set(mesh_df[mesh_col].dropna().astype(str).unique())
        print(f"  Loaded {len(valid_mesh_ids):,} valid MESH depth 3/4/5 identifiers")
        
        return valid_mesh_ids
    else:
        print(f"  ⚠️ Warning: MESH filter file not found at {mesh_file}")
        print(f"  Using dummy filter - this will affect results!")
        return set()

def filter_by_mesh_depth(df, valid_mesh_ids):
    """Filter to keep only OMIM entries OR MESH entries in depth 3/4/5 list."""
    print("Filtering by MESH depth 3/4/5 and OMIM...")
    
    initial_count = len(df)
    initial_pmids = df['pmid'].nunique()
    
    # Keep if: (1) OMIM entry, OR (2) MESH entry in our depth 3/4/5 list
    mask_omim = df['is_omim'] == True
    mask_mesh_valid = df['mesh_id_clean'].isin(valid_mesh_ids)
    
    # Combine conditions
    keep_mask = mask_omim | mask_mesh_valid
    df_filtered = df[keep_mask].copy()
    
    final_count = len(df_filtered)
    final_pmids = df_filtered['pmid'].nunique()
    
    omim_kept = df_filtered['is_omim'].sum()
    mesh_kept = final_count - omim_kept
    
    print(f"  Initial records: {initial_count:,}")
    print(f"  Filtered records: {final_count:,} ({final_count/initial_count*100:.2f}% retained)")
    print(f"  Initial PMIDs: {initial_pmids:,}")
    print(f"  Filtered PMIDs: {final_pmids:,} ({final_pmids/initial_pmids*100:.2f}% retained)")
    print(f"  OMIM entries kept: {omim_kept:,}")
    print(f"  MESH 3/4/5 entries kept: {mesh_kept:,}")
    
    return df_filtered, {
        'initial_records': initial_count,
        'final_records': final_count,
        'initial_pmids': initial_pmids,
        'final_pmids': final_pmids,
        'retention_rate': final_count/initial_count,
        'pmid_retention': final_pmids/initial_pmids,
        'omim_kept': omim_kept,
        'mesh_kept': mesh_kept
    }

def create_unique_pmid_list(df):
    """Create unique PMID list for downstream processing."""
    print("Creating unique PMID list...")
    
    # Get unique PMIDs
    unique_pmids = df[['pmid']].drop_duplicates()
    
    print(f"  Unique PMIDs: {len(unique_pmids):,}")
    
    return unique_pmids

def save_results(df_pmids, output_path, filter_stats):
    """Save results with comprehensive metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    df_pmids.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df_pmids):,} unique PMIDs")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase1_345_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 1: Disease-Relevant Publication Filter (3/4/5 Levels) Processing Statistics\\n")
        f.write("=" * 80 + "\\n\\n")
        
        f.write("Key Change: Uses MESH depth 3/4/5 levels (not just 3/4)\\n\\n")
        
        f.write("Input files:\\n")
        f.write("  - cleaned/disease_pmid_mesh.txt (disease PubTator data)\\n")
        f.write("  - rawdata/charlie/mesh_depth_3_4_5_no_explode.dta (MESH filter)\\n\\n")
        
        f.write(f"Output file: {output_path.name}\\n\\n")
        
        f.write("STEP 1: MESH/OMIM Processing\\n")
        f.write("-" * 27 + "\\n")
        f.write(f"Initial disease records: {filter_stats['initial_records']:,}\\n")
        f.write(f"OMIM entries kept: {filter_stats['omim_kept']:,}\\n")
        f.write(f"MESH 3/4/5 entries kept: {filter_stats['mesh_kept']:,}\\n")
        f.write(f"Total records retained: {filter_stats['final_records']:,}\\n")
        f.write(f"Retention rate: {filter_stats['retention_rate']*100:.2f}%\\n\\n")
        
        f.write("STEP 2: Unique PMID Extraction\\n")
        f.write("-" * 30 + "\\n")
        f.write(f"Initial PMIDs: {filter_stats['initial_pmids']:,}\\n")
        f.write(f"Disease-relevant PMIDs: {filter_stats['final_pmids']:,}\\n")
        f.write(f"PMID retention: {filter_stats['pmid_retention']*100:.2f}%\\n\\n")
        
        f.write("Final Output Summary\\n")
        f.write("-" * 20 + "\\n")
        f.write(f"Unique PMIDs: {len(df_pmids):,}\\n")
        f.write(f"Output file size: {file_size:.1f} MB\\n\\n")
        
        f.write("Next Steps:\\n")
        f.write("  1. Use these PMIDs to filter gene PubTator data\\n")
        f.write("  2. Continue with gene-disease intersection analysis\\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    disease_input = base_dir / "cleaned" / "disease_pmid_mesh.txt"
    output_file = base_dir / "processed_345" / "disease_relevant_pmids_345.parquet"
    
    print("=== Phase 1: Disease Filter (3/4/5 Levels) ===\\n")
    
    # Step 1: Load disease PubTator data
    df_disease = load_disease_data(disease_input)
    
    # Step 2: Process MESH/OMIM identifiers
    print("\\n" + "="*50)
    df_processed = process_mesh_identifiers(df_disease)
    
    # Step 3: Load MESH depth 3/4/5 filter
    print("\\n" + "="*50)
    valid_mesh_ids = load_mesh_depth_filter("mesh_depth_3_4_5")
    
    # Step 4: Filter by MESH depth and OMIM
    print("\\n" + "="*50)
    df_filtered, filter_stats = filter_by_mesh_depth(df_processed, valid_mesh_ids)
    
    # Step 5: Create unique PMID list
    print("\\n" + "="*50)
    df_pmids = create_unique_pmid_list(df_filtered)
    
    # Step 6: Save results
    print("\\n" + "="*50)
    save_results(df_pmids, output_file, filter_stats)
    
    print(f"\\n=== Phase 1 (3/4/5) Complete ===")
    print(f"Output: {output_file}")
    print(f"Disease-relevant PMIDs (3/4/5 levels): {len(df_pmids):,}")

if __name__ == "__main__":
    main()