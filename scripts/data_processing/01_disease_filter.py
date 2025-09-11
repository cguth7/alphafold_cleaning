#!/usr/bin/env python3
"""
Phase 1: Disease-Relevant Publication Filter

Filter disease PubTator data to keep only PMIDs with:
1. MESH diseases at depth levels 3/4 OR
2. OMIM diseases

Input:  cleaned/disease_pmid_mesh.txt (163M rows, 3.4GB)  
Output: processed/disease_relevant_pmids.parquet

Processing steps:
1. Load cleaned disease data (pmid, mesh_id)
2. Strip MESH:/OMIM: prefixes and flag OMIM entries
3. Filter by MESH depth 3/4 reference list  
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
    
    # Verify cleaning worked
    remaining_prefixes = df['mesh_id_clean'].str.contains('MESH:|OMIM:', na=False).sum()
    if remaining_prefixes > 0:
        print(f"  WARNING: {remaining_prefixes:,} entries still have prefixes!")
    
    print(f"  Cleaned {len(df):,} mesh identifiers")
    
    return df

def load_mesh_depth_filter(file_path):
    """Load MESH depth 3/4 reference data."""
    print(f"Loading MESH depth filter from: {file_path}")
    
    df = pd.read_stata(file_path)
    
    print(f"  Loaded {len(df):,} MESH descriptors")
    print(f"  Depth distribution:")
    print(df['depth'].value_counts().sort_index().to_string())
    
    # Keep only unique MESH IDs (remove depth duplicates)
    mesh_ids = df['DescriptorUI'].unique()
    print(f"  Unique MESH IDs: {len(mesh_ids):,}")
    
    return set(mesh_ids)

def filter_by_mesh_depth(df, mesh_depth_set):
    """Filter disease data by MESH depth 3/4 or OMIM."""
    print("Filtering by MESH depth 3/4 or OMIM...")
    
    initial_count = len(df)
    
    # Check which entries match MESH depth filter
    df['in_mesh_depth'] = df['mesh_id_clean'].isin(mesh_depth_set)
    mesh_matches = df['in_mesh_depth'].sum()
    
    # Keep either OMIM entries OR MESH depth matches
    df_filtered = df[df['is_omim'] | df['in_mesh_depth']].copy()
    
    final_count = len(df_filtered)
    omim_kept = df_filtered['is_omim'].sum()
    mesh_kept = df_filtered['in_mesh_depth'].sum()
    
    print(f"  Initial rows: {initial_count:,}")
    print(f"  MESH depth matches: {mesh_matches:,} ({mesh_matches/initial_count*100:.2f}%)")
    print(f"  Kept - OMIM: {omim_kept:,}, MESH depth: {mesh_kept:,}")
    print(f"  Final rows: {final_count:,} ({final_count/initial_count*100:.2f}% retention)")
    
    # Sample of excluded entries for documentation
    excluded_sample = df[~(df['is_omim'] | df['in_mesh_depth'])]['mesh_id'].head(5).tolist()
    print(f"  Sample excluded: {excluded_sample}")
    
    return df_filtered, {
        'initial_count': initial_count,
        'final_count': final_count, 
        'omim_kept': omim_kept,
        'mesh_kept': mesh_kept,
        'retention_rate': final_count/initial_count,
        'excluded_sample': excluded_sample
    }

def create_unique_pmid_list(df):
    """Create unique PMID list from filtered data."""
    print("Creating unique PMID list...")
    
    initial_rows = len(df)
    unique_pmids = df['pmid'].drop_duplicates().sort_values()
    unique_count = len(unique_pmids)
    
    print(f"  Initial disease mentions: {initial_rows:,}")
    print(f"  Unique PMIDs: {unique_count:,}")
    print(f"  Avg mentions per paper: {initial_rows/unique_count:.2f}")
    
    return unique_pmids.to_frame()

def save_results(pmid_df, output_path, stats):
    """Save results with metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    pmid_df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(pmid_df):,} unique PMIDs")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase1_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 1: Disease Filter Processing Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input file: cleaned/disease_pmid_mesh.txt\n")
        f.write(f"Output file: {output_path.name}\n\n")
        f.write(f"Initial disease mentions: {stats['initial_count']:,}\n")
        f.write(f"Final disease mentions: {stats['final_count']:,}\n")
        f.write(f"Retention rate: {stats['retention_rate']*100:.2f}%\n\n")
        f.write(f"Kept entries:\n")
        f.write(f"  OMIM diseases: {stats['omim_kept']:,}\n") 
        f.write(f"  MESH depth 3/4: {stats['mesh_kept']:,}\n\n")
        f.write(f"Unique PMIDs: {len(pmid_df):,}\n")
        f.write(f"Output file size: {file_size:.1f} MB\n\n")
        f.write(f"Sample excluded identifiers:\n")
        for i, ex in enumerate(stats['excluded_sample'], 1):
            f.write(f"  {i}. {ex}\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    disease_input = base_dir / "cleaned" / "disease_pmid_mesh.txt"
    mesh_filter = base_dir / "raw_data" / "mesh_levels" / "mesh_depth_3_4_no_explode.dta" 
    output_file = base_dir / "processed" / "disease_relevant_pmids.parquet"
    
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print("=== Phase 1: Disease-Relevant Publication Filter ===\n")
    
    # Step 1: Load disease data
    df_disease = load_disease_data(disease_input)
    
    # Step 2: Process identifiers
    df_disease = process_mesh_identifiers(df_disease)
    
    # Step 3: Load MESH depth filter
    mesh_depth_set = load_mesh_depth_filter(mesh_filter)
    
    # Step 4: Filter by MESH depth or OMIM
    df_filtered, stats = filter_by_mesh_depth(df_disease, mesh_depth_set)
    
    # Step 5: Create unique PMID list
    pmid_df = create_unique_pmid_list(df_filtered)
    
    # Step 6: Save results
    save_results(pmid_df, output_file, stats)
    
    print(f"\n=== Phase 1 Complete ===")
    print(f"Output: {output_file}")
    print(f"Unique disease-relevant PMIDs: {len(pmid_df):,}")

if __name__ == "__main__":
    main()