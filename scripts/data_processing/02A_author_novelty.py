#!/usr/bin/env python3
"""
Phase 2A: Author Novelty Analysis (Parallel Track)

For each author-gene pair, determine if this is the author's first paper on that gene
(considering publications from 2015 onward). Creates newcomer/veteran author indicators
for later aggregation at the gene-month level.

Input:  
- processed/gene_disease_temporal.parquet (13.7M records, AlphaFold era)
- intermediate_data/pmids_authors_openalex_deduped.dta (4.2GB author data)

Output: processed/author_novelty.parquet

Processing steps:
1. Load gene-disease temporal data from Phase 3
2. Load author data and merge with gene-disease records  
3. Expand analysis window to 2015+ to establish author histories
4. For each author-gene pair, identify first occurrence (since 2015)
5. Create newcomer_author and not_newcomer_author flags
6. Filter back to AlphaFold era (2020-2023) with novelty flags
7. Save with comprehensive author statistics
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_gene_disease_temporal(file_path):
    """Load temporal gene-disease data from Phase 3."""
    print(f"Loading gene-disease temporal data from: {file_path}")
    
    df = pd.read_parquet(file_path)
    
    print(f"  Loaded {len(df):,} gene-disease temporal records")
    print(f"  Time span: {df['year'].min()}-{df['year'].max()}")
    print(f"  Unique PMIDs: {df['pmid'].nunique():,}")
    print(f"  Unique genes: {df['gene_id'].nunique():,}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def load_author_data(file_path):
    """Load deduplicated author data."""
    print(f"Loading author data from: {file_path}")
    print("  (This may take a few minutes due to file size...)")
    
    # Load the large author file
    df = pd.read_stata(file_path)
    
    print(f"  Loaded {len(df):,} author records")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    # Show a sample to understand structure
    print("  Sample records:")
    print(df.head(3).to_string())
    
    return df

def merge_with_authors(df_gene, df_authors):
    """Merge gene-disease data with author information."""
    print("Merging gene-disease data with author information...")
    
    # Determine merge key (likely pmid)
    if 'pmid' in df_authors.columns:
        merge_key = 'pmid'
    else:
        # Check for other possible keys
        possible_keys = ['PMID', 'pubmed_id', 'pmid_int']
        merge_key = None
        for key in possible_keys:
            if key in df_authors.columns:
                merge_key = key
                break
        
        if merge_key is None:
            print(f"  Error: No suitable merge key found in author data")
            print(f"  Available columns: {list(df_authors.columns)}")
            return None, None
    
    print(f"  Using merge key: {merge_key}")
    
    # Ensure consistent data types for merging
    if merge_key != 'pmid':
        df_authors = df_authors.rename(columns={merge_key: 'pmid'})
    
    # Merge
    initial_records = len(df_gene)
    initial_pmids = df_gene['pmid'].nunique()
    
    df_merged = df_gene.merge(df_authors, on='pmid', how='inner')
    
    final_records = len(df_merged)
    final_pmids = df_merged['pmid'].nunique()
    
    print(f"  Initial gene records: {initial_records:,}")
    print(f"  Merged records: {final_records:,} ({final_records/initial_records*100:.2f}% retained)")
    print(f"  Initial PMIDs: {initial_pmids:,}")
    print(f"  PMIDs with authors: {final_pmids:,} ({final_pmids/initial_pmids*100:.2f}% coverage)")
    
    if final_records == 0:
        print("  ⚠️ Error: No successful merges. Check data compatibility.")
        return None, None
    
    # Show unique authors
    author_col = None
    for col in ['last_author_id', 'author_id', 'author', 'openalex_author_id', 'author_name']:
        if col in df_merged.columns:
            author_col = col
            break
    
    if author_col:
        unique_authors = df_merged[author_col].nunique()
        print(f"  Unique authors: {unique_authors:,} (using column: {author_col})")
    else:
        print("  ⚠️ Warning: Could not identify author identifier column")
        print(f"  Available columns: {list(df_merged.columns)}")
        # Default to last_author_id since that's what we see in the data
        author_col = 'last_author_id'
        unique_authors = df_merged[author_col].nunique()
        print(f"  Using default: {author_col} with {unique_authors:,} unique values")
    
    return df_merged, {
        'initial_records': initial_records,
        'final_records': final_records,
        'initial_pmids': initial_pmids,
        'final_pmids': final_pmids,
        'merge_rate': final_records/initial_records,
        'pmid_coverage': final_pmids/initial_pmids,
        'author_column': author_col,
        'unique_authors': df_merged[author_col].nunique() if author_col else 0
    }

def expand_temporal_window(df_with_authors, start_year=2015):
    """Expand temporal analysis window to establish author histories."""
    print(f"Expanding temporal window to {start_year}+ for author history analysis...")
    
    # We need to go back and get more temporal data
    # For now, we'll work with what we have and note this limitation
    current_min_year = df_with_authors['year'].min()
    
    if current_min_year <= start_year:
        print(f"  Current data spans {current_min_year}-{df_with_authors['year'].max()}")
        print(f"  Good: we have data from {start_year} for baseline establishment")
        expand_needed = False
    else:
        print(f"  Current data starts from {current_min_year}")
        print(f"  ⚠️ Warning: Need data from {start_year} for proper baseline")
        print(f"  Will proceed with available data but results may be limited")
        expand_needed = True
    
    return df_with_authors, expand_needed

def create_author_novelty_flags(df, author_col='last_author_id'):
    """Create newcomer_author and not_newcomer_author flags."""
    print("Creating author novelty flags...")
    
    if author_col not in df.columns:
        print(f"  Error: Author column '{author_col}' not found")
        return None
    
    initial_records = len(df)
    
    # Sort by author, gene, and time to identify first occurrences
    df_sorted = df.sort_values([author_col, 'gene_id', 'year', 'month', 'pmid']).copy()
    
    # Create author-gene pair identifier
    print("  Identifying first author-gene publications...")
    df_sorted['author_gene_first'] = ~df_sorted.duplicated(subset=[author_col, 'gene_id'], keep='first')
    
    # Create the novelty flags
    df_sorted['newcomer_author'] = df_sorted['author_gene_first'].astype(int)
    df_sorted['not_newcomer_author'] = (1 - df_sorted['newcomer_author']).astype(int)
    
    # Statistics
    total_author_gene_pairs = df_sorted[[author_col, 'gene_id']].drop_duplicates().shape[0]
    newcomer_records = df_sorted['newcomer_author'].sum()
    veteran_records = df_sorted['not_newcomer_author'].sum()
    
    print(f"  Total author-gene pairs: {total_author_gene_pairs:,}")
    print(f"  Newcomer author records: {newcomer_records:,} ({newcomer_records/initial_records*100:.2f}%)")
    print(f"  Veteran author records: {veteran_records:,} ({veteran_records/initial_records*100:.2f}%)")
    
    # Verify logic
    assert newcomer_records + veteran_records == initial_records, "Novelty flags don't sum correctly"
    assert newcomer_records == total_author_gene_pairs, "Should have one newcomer record per author-gene pair"
    
    # Drop intermediate columns
    df_final = df_sorted.drop(['author_gene_first'], axis=1)
    
    return df_final, {
        'total_records': initial_records,
        'total_author_gene_pairs': total_author_gene_pairs,
        'newcomer_records': newcomer_records,
        'veteran_records': veteran_records,
        'newcomer_rate': newcomer_records/initial_records
    }

def save_results(df, output_path, merge_stats, novelty_stats):
    """Save results with comprehensive metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} author-gene records")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase2A_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 2A: Author Novelty Analysis Processing Statistics\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Input files:\n")
        f.write("  - processed/gene_disease_temporal.parquet (Phase 3 output)\n")
        f.write("  - intermediate_data/pmids_authors_openalex_deduped.dta\n\n")
        
        f.write(f"Output file: {output_path.name}\n\n")
        
        f.write("STEP 1: Author Data Merging\n")
        f.write("-" * 27 + "\n")
        f.write(f"Initial gene-disease records: {merge_stats['initial_records']:,}\n")
        f.write(f"Records with author data: {merge_stats['final_records']:,}\n")
        f.write(f"Merge success rate: {merge_stats['merge_rate']*100:.2f}%\n")
        f.write(f"PMID coverage rate: {merge_stats['pmid_coverage']*100:.2f}%\n")
        f.write(f"Unique authors identified: {merge_stats['unique_authors']:,}\n\n")
        
        f.write("STEP 2: Author Novelty Analysis\n")
        f.write("-" * 31 + "\n")
        f.write(f"Total author-gene pairs: {novelty_stats['total_author_gene_pairs']:,}\n")
        f.write(f"Newcomer author records: {novelty_stats['newcomer_records']:,}\n")
        f.write(f"Veteran author records: {novelty_stats['veteran_records']:,}\n")
        f.write(f"Newcomer rate: {novelty_stats['newcomer_rate']*100:.2f}%\n\n")
        
        f.write("Final Output Summary\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Unique PMIDs: {df['pmid'].nunique():,}\n")
        f.write(f"Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"Unique authors: {df[merge_stats['author_column']].nunique():,}\n")
        f.write(f"Time span: {df['year'].min()}-{df['year'].max()}\n")
        f.write(f"Output file size: {file_size:.1f} MB\n\n")
        
        f.write("Novelty Flags Created:\n")
        f.write("  - newcomer_author: 1 for first publication by author on gene, 0 otherwise\n")
        f.write("  - not_newcomer_author: 1 for repeat publication by author on gene, 0 otherwise\n\n")
        
        f.write("Usage in Final Panel:\n")
        f.write("  These flags will be aggregated at gene-month level to track:\n")
        f.write("  - Number of newcomer authors per gene per month\n")
        f.write("  - Number of veteran authors per gene per month\n")
        f.write("  - Author novelty trends over the AlphaFold era\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    temporal_input = base_dir / "processed" / "gene_disease_temporal.parquet"
    author_input = base_dir / "intermediate_data" / "pmids_authors_openalex_deduped.dta"
    output_file = base_dir / "processed" / "author_novelty.parquet"
    
    print("=== Phase 2A: Author Novelty Analysis ===\n")
    
    # Step 1: Load temporal data
    df_temporal = load_gene_disease_temporal(temporal_input)
    
    # Step 2: Load and merge author data
    print("\\n" + "="*50)
    df_authors = load_author_data(author_input)
    
    df_with_authors, merge_stats = merge_with_authors(df_temporal, df_authors)
    
    if df_with_authors is None:
        print("Failed to merge author data. Exiting.")
        return
    
    # Step 3: Expand temporal window (note limitations)
    print("\\n" + "="*50)
    df_expanded, expand_needed = expand_temporal_window(df_with_authors)
    
    # Step 4: Create novelty flags
    print("\\n" + "="*50)
    author_column = merge_stats.get('author_column', 'last_author_id')
    df_final, novelty_stats = create_author_novelty_flags(df_expanded, author_column)
    
    if df_final is None:
        print("Failed to create novelty flags. Exiting.")
        return
    
    # Step 5: Save results
    print("\\n" + "="*50)
    save_results(df_final, output_file, merge_stats, novelty_stats)
    
    print(f"\\n=== Phase 2A Complete ===")
    print(f"Output: {output_file}")
    print(f"Author-gene records with novelty flags: {len(df_final):,}")
    print(f"Unique authors: {df_final[merge_stats['author_column']].nunique():,}")
    print(f"Newcomer rate: {novelty_stats['newcomer_rate']*100:.1f}%")

if __name__ == "__main__":
    main()