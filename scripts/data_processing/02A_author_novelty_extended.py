#!/usr/bin/env python3
"""
Phase 2A Extended: Author Novelty Analysis with Proper Baseline

For each author-gene pair, determine if this is the author's first paper on that gene
(considering publications from 2015 onward). Creates newcomer/veteran author indicators
for later aggregation at the gene-month level, but uses proper 2015+ baseline.

Input:  
- processed/gene_disease_temporal_extended.parquet (22.2M records, extended period)
- intermediate_data/pmids_authors_openalex_deduped.dta (4.2GB author data)

Output: processed/author_novelty_extended.parquet

Processing steps:
1. Load extended gene-disease temporal data (2015-2023)
2. Load author data and merge with gene-disease records  
3. For each author-gene pair, identify first occurrence (since 2015)
4. Create newcomer_author and not_newcomer_author flags
5. Filter back to AlphaFold era (2020-2023) with proper novelty flags
6. Save with comprehensive author statistics
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_extended_temporal_data(file_path):
    """Load extended gene-disease temporal data."""
    print(f"Loading extended temporal data from: {file_path}")
    
    df = pd.read_parquet(file_path)
    
    print(f"  Loaded {len(df):,} extended temporal records")
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
    print("Merging extended temporal data with author information...")
    
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

def create_author_novelty_flags_extended(df, author_col='last_author_id'):
    """Create newcomer_author and not_newcomer_author flags using full 2015+ baseline."""
    print("Creating author novelty flags with 2015+ baseline...")
    
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
    
    # Statistics across full period
    total_author_gene_pairs = df_sorted[[author_col, 'gene_id']].drop_duplicates().shape[0]
    newcomer_records = df_sorted['newcomer_author'].sum()
    veteran_records = df_sorted['not_newcomer_author'].sum()
    
    print(f"  Total author-gene pairs (2015+): {total_author_gene_pairs:,}")
    print(f"  Newcomer author records: {newcomer_records:,} ({newcomer_records/initial_records*100:.2f}%)")
    print(f"  Veteran author records: {veteran_records:,} ({veteran_records/initial_records*100:.2f}%)")
    
    # Show breakdown by year
    yearly_stats = df_sorted.groupby('year').agg({
        'newcomer_author': 'sum',
        'not_newcomer_author': 'sum'
    })
    print(f"  Yearly breakdown:")
    for year, row in yearly_stats.iterrows():
        total_year = row['newcomer_author'] + row['not_newcomer_author'] 
        if total_year > 0:
            newcomer_pct = row['newcomer_author'] / total_year * 100
            print(f"    {year}: {row['newcomer_author']:,} newcomers ({newcomer_pct:.1f}%), {row['not_newcomer_author']:,} veterans")
    
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
        'newcomer_rate': newcomer_records/initial_records,
        'yearly_stats': yearly_stats
    }

def filter_to_alphafold_era(df):
    """Filter to AlphaFold era while preserving novelty flags established from 2015+ baseline."""
    print("Filtering to AlphaFold era (2020-2023) while preserving extended baseline...")
    
    initial_records = len(df)
    initial_years = df['year'].value_counts().sort_index()
    
    print(f"  Pre-filter breakdown:")
    for year, count in initial_years.items():
        print(f"    {year}: {count:,} records")
    
    # Filter to 2020-2023
    df_alphafold = df[df['year'].between(2020, 2023)].copy()
    
    final_records = len(df_alphafold)
    final_years = df_alphafold['year'].value_counts().sort_index()
    
    print(f"  Post-filter (AlphaFold era):")
    print(f"    Total records: {final_records:,} ({final_records/initial_records*100:.2f}% of extended data)")
    for year, count in final_years.items():
        print(f"    {year}: {count:,} records")
    
    # Show novelty stats in AlphaFold era
    alphafold_newcomer = df_alphafold['newcomer_author'].sum()
    alphafold_veteran = df_alphafold['not_newcomer_author'].sum()
    
    print(f"  AlphaFold era novelty (with 2015+ baseline):")
    print(f"    Newcomer papers: {alphafold_newcomer:,} ({alphafold_newcomer/final_records*100:.1f}%)")
    print(f"    Veteran papers: {alphafold_veteran:,} ({alphafold_veteran/final_records*100:.1f}%)")
    
    return df_alphafold, {
        'initial_records': initial_records,
        'final_records': final_records,
        'retention_rate': final_records/initial_records,
        'alphafold_newcomer_rate': alphafold_newcomer/final_records if final_records > 0 else 0
    }

def save_results(df, output_path, merge_stats, novelty_stats, filter_stats):
    """Save results with comprehensive metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} author-gene records with extended baseline")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase2A_extended_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 2A Extended: Author Novelty Analysis Processing Statistics\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("Input files:\n")
        f.write("  - processed/gene_disease_temporal_extended.parquet (2015-2023 data)\n")
        f.write("  - intermediate_data/pmids_authors_openalex_deduped.dta\n\n")
        
        f.write(f"Output file: {output_path.name}\n\n")
        
        f.write("STEP 1: Author Data Merging\n")
        f.write("-" * 27 + "\n")
        f.write(f"Initial extended temporal records: {merge_stats['initial_records']:,}\n")
        f.write(f"Records with author data: {merge_stats['final_records']:,}\n")
        f.write(f"Merge success rate: {merge_stats['merge_rate']*100:.2f}%\n")
        f.write(f"PMID coverage rate: {merge_stats['pmid_coverage']*100:.2f}%\n")
        f.write(f"Unique authors identified: {merge_stats['unique_authors']:,}\n\n")
        
        f.write("STEP 2: Author Novelty Analysis (2015+ Baseline)\n")
        f.write("-" * 46 + "\n")
        f.write(f"Total author-gene pairs: {novelty_stats['total_author_gene_pairs']:,}\n")
        f.write(f"Newcomer author records: {novelty_stats['newcomer_records']:,}\n")
        f.write(f"Veteran author records: {novelty_stats['veteran_records']:,}\n")
        f.write(f"Newcomer rate: {novelty_stats['newcomer_rate']*100:.2f}%\n\n")
        
        f.write("STEP 3: AlphaFold Era Filtering\n")
        f.write("-" * 30 + "\n")
        f.write(f"Extended baseline records: {filter_stats['initial_records']:,}\n")
        f.write(f"AlphaFold era records: {filter_stats['final_records']:,}\n")
        f.write(f"Retention rate: {filter_stats['retention_rate']*100:.2f}%\n")
        f.write(f"AlphaFold newcomer rate: {filter_stats['alphafold_newcomer_rate']*100:.2f}%\n\n")
        
        f.write("Final Output Summary\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Unique PMIDs: {df['pmid'].nunique():,}\n")
        f.write(f"Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"Unique authors: {df[merge_stats['author_column']].nunique():,}\n")
        f.write(f"Time span: {df['year'].min()}-{df['year'].max()}\n")
        f.write(f"Output file size: {file_size:.1f} MB\n\n")
        
        f.write("Key Improvement:\n")
        f.write("  - Uses proper 2015+ baseline for author novelty classification\n")
        f.write("  - More accurate newcomer/veteran distinctions in AlphaFold era\n")
        f.write("  - Addresses the original issue in 02A analysis\n\n")
        
        f.write("Novelty Flags Created:\n")
        f.write("  - newcomer_author: 1 for first publication by author on gene since 2015, 0 otherwise\n")
        f.write("  - not_newcomer_author: 1 for repeat publication by author on gene, 0 otherwise\n\n")
        
        f.write("Next Steps:\n")
        f.write("  - Aggregate to gene-month level for panel analysis\n")
        f.write("  - Merge aggregated author metrics into final panel\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    temporal_input = base_dir / "processed" / "gene_disease_temporal_extended.parquet"
    author_input = base_dir / "intermediate_data" / "pmids_authors_openalex_deduped.dta"
    output_file = base_dir / "processed" / "author_novelty_extended.parquet"
    
    print("=== Phase 2A Extended: Author Novelty Analysis ===\n")
    
    # Step 1: Load extended temporal data
    df_temporal = load_extended_temporal_data(temporal_input)
    
    # Step 2: Load and merge author data
    print("\\n" + "="*50)
    df_authors = load_author_data(author_input)
    
    df_with_authors, merge_stats = merge_with_authors(df_temporal, df_authors)
    
    if df_with_authors is None:
        print("Failed to merge author data. Exiting.")
        return
    
    # Step 3: Create novelty flags with extended baseline
    print("\\n" + "="*50)
    author_column = merge_stats.get('author_column', 'last_author_id')
    df_with_novelty, novelty_stats = create_author_novelty_flags_extended(df_with_authors, author_column)
    
    if df_with_novelty is None:
        print("Failed to create novelty flags. Exiting.")
        return
    
    # Step 4: Filter to AlphaFold era while preserving extended baseline
    print("\\n" + "="*50)
    df_final, filter_stats = filter_to_alphafold_era(df_with_novelty)
    
    # Step 5: Save results
    print("\\n" + "="*50)
    save_results(df_final, output_file, merge_stats, novelty_stats, filter_stats)
    
    print(f"\\n=== Phase 2A Extended Complete ===")
    print(f"Output: {output_file}")
    print(f"AlphaFold era records with proper 2015+ baseline: {len(df_final):,}")
    print(f"Unique authors: {df_final[merge_stats['author_column']].nunique():,}")
    print(f"Improved newcomer rate: {filter_stats['alphafold_newcomer_rate']*100:.1f}%")

if __name__ == "__main__":
    main()