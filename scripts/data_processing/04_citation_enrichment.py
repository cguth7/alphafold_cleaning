#!/usr/bin/env python3
"""
Phase 4: Citation Enrichment

Add citation data and create quality percentile flags at multiple temporal granularities:
- Annual percentiles (top 25%, 10%, 5% within each year)
- Quarterly percentiles (top 25%, 10%, 5% within each quarter)  
- Bi-monthly percentiles (top 25%, 10%, 5% within each bi-month)

Input:  
- processed/author_novelty.parquet (13.7M records with author info)
- raw_data/icite_cites_join.dta (6.7M PMIDs with citation counts)

Output: processed/citation_enriched.parquet

Processing steps:
1. Load author novelty data from Phase 2A
2. Merge with iCite citation data
3. Create paper-level aggregations with max citation flags per PMID
4. Calculate percentiles within year/quarter/bi-month periods
5. Create high-citation flags at each temporal granularity  
6. Merge flags back to gene-level records
7. Save with comprehensive citation statistics
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_author_novelty_data(file_path):
    """Load author novelty data from Phase 2A."""
    print(f"Loading author novelty data from: {file_path}")
    
    df = pd.read_parquet(file_path)
    
    print(f"  Loaded {len(df):,} author-gene records")
    print(f"  Time span: {df['year'].min()}-{df['year'].max()}")
    print(f"  Unique PMIDs: {df['pmid'].nunique():,}")
    print(f"  Unique genes: {df['gene_id'].nunique():,}")
    print(f"  Unique authors: {df['last_author_id'].nunique():,}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def load_citation_data(file_path):
    """Load iCite citation data."""
    print(f"Loading citation data from: {file_path}")
    
    df = pd.read_stata(file_path)
    
    print(f"  Loaded {len(df):,} PMID-citation records")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Citation range: {df['icite_cites'].min()}-{df['icite_cites'].max()}")
    print(f"  Mean citations: {df['icite_cites'].mean():.2f}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    # Data quality checks
    null_cites = df['icite_cites'].isnull().sum()
    if null_cites > 0:
        print(f"  Found {null_cites:,} records with missing citations")
        df['icite_cites'] = df['icite_cites'].fillna(0)
        print(f"  Filled missing citations with 0")
    
    return df

def merge_with_citations(df_author, df_citations):
    """Merge author data with citation information."""
    print("Merging author data with citation information...")
    
    initial_records = len(df_author)
    initial_pmids = df_author['pmid'].nunique()
    
    # Merge on pmid - use left join to keep all author records
    df_merged = df_author.merge(df_citations, on='pmid', how='left')
    
    final_records = len(df_merged)
    final_pmids = df_merged['pmid'].nunique()
    
    # Fill missing citations with 0 (papers not in iCite database)
    df_merged['icite_cites'] = df_merged['icite_cites'].fillna(0)
    
    # Statistics
    pmids_with_cites = (df_merged['icite_cites'] > 0).groupby(df_merged['pmid']).any().sum()
    
    print(f"  Initial records: {initial_records:,}")
    print(f"  Final records: {final_records:,} (should be same - left join)")
    print(f"  PMIDs with citation data: {pmids_with_cites:,}/{final_pmids:,} ({pmids_with_cites/final_pmids*100:.2f}%)")
    
    cite_stats = df_merged['icite_cites'].describe()
    print(f"  Citation statistics:")
    print(f"    Mean: {cite_stats['mean']:.2f}")
    print(f"    Median: {cite_stats['50%']:.2f}")
    print(f"    Max: {cite_stats['max']:.0f}")
    
    return df_merged, {
        'initial_records': initial_records,
        'final_records': final_records,
        'pmids_with_citations': pmids_with_cites,
        'total_pmids': final_pmids,
        'citation_coverage': pmids_with_cites/final_pmids,
        'mean_citations': cite_stats['mean'],
        'median_citations': cite_stats['50%']
    }

def create_paper_level_data(df):
    """Create paper-level data with max citation flags per PMID."""
    print("Creating paper-level citation aggregations...")
    
    # Group by PMID and take the first occurrence for metadata, max for any gene-level flags
    paper_data = (df.groupby('pmid')
                    .agg({
                        'year': 'first',
                        'month': 'first', 
                        'ym': 'first',
                        'quarter': 'first',
                        'year_quarter': 'first',
                        'bimonth': 'first',
                        'year_bimonth': 'first',
                        'icite_cites': 'first',  # Should be same for all genes in same paper
                        'newcomer_author': 'max',  # Any newcomer authors in this paper
                        'not_newcomer_author': 'max'  # Any veteran authors in this paper
                    })
                    .reset_index())
    
    print(f"  Created {len(paper_data):,} paper-level records")
    print(f"  Papers with newcomer authors: {paper_data['newcomer_author'].sum():,} ({paper_data['newcomer_author'].mean()*100:.1f}%)")
    print(f"  Papers with veteran authors: {paper_data['not_newcomer_author'].sum():,} ({paper_data['not_newcomer_author'].mean()*100:.1f}%)")
    
    return paper_data

def calculate_citation_percentiles(paper_data):
    """Calculate citation percentiles within year/quarter/bi-month periods."""
    print("Calculating citation percentiles...")
    
    # Annual percentiles
    print("  Calculating annual percentiles...")
    paper_data = paper_data.merge(
        paper_data.groupby('year')['icite_cites'].quantile([0.75, 0.90, 0.95]).unstack().rename(columns={
            0.75: 'p75_year', 0.90: 'p90_year', 0.95: 'p95_year'
        }),
        on='year'
    )
    
    paper_data['hi25_year'] = (paper_data['icite_cites'] >= paper_data['p75_year']).astype(int)
    paper_data['hi10_year'] = (paper_data['icite_cites'] >= paper_data['p90_year']).astype(int)  
    paper_data['hi05_year'] = (paper_data['icite_cites'] >= paper_data['p95_year']).astype(int)
    
    # Quarterly percentiles
    print("  Calculating quarterly percentiles...")
    paper_data = paper_data.merge(
        paper_data.groupby('year_quarter')['icite_cites'].quantile([0.75, 0.90, 0.95]).unstack().rename(columns={
            0.75: 'p75_quarter', 0.90: 'p90_quarter', 0.95: 'p95_quarter'
        }),
        on='year_quarter'
    )
    
    paper_data['hi25_quarter'] = (paper_data['icite_cites'] >= paper_data['p75_quarter']).astype(int)
    paper_data['hi10_quarter'] = (paper_data['icite_cites'] >= paper_data['p90_quarter']).astype(int)
    paper_data['hi05_quarter'] = (paper_data['icite_cites'] >= paper_data['p95_quarter']).astype(int)
    
    # Bi-monthly percentiles  
    print("  Calculating bi-monthly percentiles...")
    paper_data = paper_data.merge(
        paper_data.groupby('year_bimonth')['icite_cites'].quantile([0.75, 0.90, 0.95]).unstack().rename(columns={
            0.75: 'p75_bimonth', 0.90: 'p90_bimonth', 0.95: 'p95_bimonth'
        }),
        on='year_bimonth'
    )
    
    paper_data['hi25_bimonth'] = (paper_data['icite_cites'] >= paper_data['p75_bimonth']).astype(int)
    paper_data['hi10_bimonth'] = (paper_data['icite_cites'] >= paper_data['p90_bimonth']).astype(int)
    paper_data['hi05_bimonth'] = (paper_data['icite_cites'] >= paper_data['p95_bimonth']).astype(int)
    
    # Summary statistics
    for period, prefix in [('year', 'year'), ('quarter', 'quarter'), ('bimonth', 'bimonth')]:
        hi25_count = paper_data[f'hi25_{prefix}'].sum()
        hi10_count = paper_data[f'hi10_{prefix}'].sum() 
        hi05_count = paper_data[f'hi05_{prefix}'].sum()
        
        print(f"  {period.title()} high-citation papers:")
        print(f"    Top 25%: {hi25_count:,} papers ({hi25_count/len(paper_data)*100:.1f}%)")
        print(f"    Top 10%: {hi10_count:,} papers ({hi10_count/len(paper_data)*100:.1f}%)")
        print(f"    Top 5%:  {hi05_count:,} papers ({hi05_count/len(paper_data)*100:.1f}%)")
    
    # Drop intermediate percentile columns
    cols_to_drop = [col for col in paper_data.columns if col.startswith(('p75_', 'p90_', 'p95_'))]
    paper_data = paper_data.drop(columns=cols_to_drop)
    
    return paper_data

def merge_citation_flags_back(df_original, paper_data):
    """Merge citation flags back to original gene-level data."""
    print("Merging citation flags back to gene-level data...")
    
    # Select citation flag columns to merge back
    citation_cols = ['pmid'] + [col for col in paper_data.columns if col.startswith('hi')]
    paper_flags = paper_data[citation_cols]
    
    # Merge back to original data
    df_enriched = df_original.merge(paper_flags, on='pmid', how='left')
    
    print(f"  Original records: {len(df_original):,}")
    print(f"  Enriched records: {len(df_enriched):,}")
    print(f"  Citation flags added: {len(citation_cols)-1}")
    
    return df_enriched

def save_results(df, output_path, merge_stats):
    """Save results with comprehensive metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} citation-enriched records")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase4_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 4: Citation Enrichment Processing Statistics\n")
        f.write("=" * 55 + "\n\n")
        
        f.write("Input files:\n")
        f.write("  - processed/author_novelty.parquet (Phase 2A output)\n")
        f.write("  - raw_data/icite_cites_join.dta (iCite citation data)\n\n")
        
        f.write(f"Output file: {output_path.name}\n\n")
        
        f.write("STEP 1: Citation Data Merging\n")
        f.write("-" * 29 + "\n")
        f.write(f"Author-novelty records: {merge_stats['initial_records']:,}\n")
        f.write(f"Records after citation merge: {merge_stats['final_records']:,}\n")
        f.write(f"PMIDs with citations: {merge_stats['pmids_with_citations']:,}/{merge_stats['total_pmids']:,}\n")
        f.write(f"Citation coverage: {merge_stats['citation_coverage']*100:.2f}%\n")
        f.write(f"Mean citations per paper: {merge_stats['mean_citations']:.2f}\n")
        f.write(f"Median citations per paper: {merge_stats['median_citations']:.2f}\n\n")
        
        f.write("STEP 2: Citation Quality Flags\n")
        f.write("-" * 30 + "\n")
        f.write("Created high-citation flags at three temporal granularities:\n")
        f.write("  - Annual: hi25_year, hi10_year, hi05_year\n")
        f.write("  - Quarterly: hi25_quarter, hi10_quarter, hi05_quarter\n") 
        f.write("  - Bi-monthly: hi25_bimonth, hi10_bimonth, hi05_bimonth\n\n")
        
        f.write("Final Output Summary\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Unique PMIDs: {df['pmid'].nunique():,}\n")
        f.write(f"Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"Unique authors: {df['last_author_id'].nunique():,}\n")
        f.write(f"Time span: {df['year'].min()}-{df['year'].max()}\n")
        f.write(f"Output file size: {file_size:.1f} MB\n\n")
        
        f.write("Ready for Phase 5: Panel Construction\n")
        f.write("  - Gene-level records with citation quality flags\n")
        f.write("  - Author novelty indicators (newcomer/veteran)\n")
        f.write("  - Multi-granularity temporal variables\n")
        f.write("  - Comprehensive metadata for aggregation\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    author_input = base_dir / "processed" / "author_novelty.parquet"
    citation_input = base_dir / "raw_data" / "icite_cites_join.dta"
    output_file = base_dir / "processed" / "citation_enriched.parquet"
    
    print("=== Phase 4: Citation Enrichment ===\n")
    
    # Step 1: Load author novelty data
    df_author = load_author_novelty_data(author_input)
    
    # Step 2: Load and merge citation data
    print("\\n" + "="*50)
    df_citations = load_citation_data(citation_input)
    
    df_with_cites, merge_stats = merge_with_citations(df_author, df_citations)
    
    # Step 3: Create paper-level aggregations
    print("\\n" + "="*50)
    paper_data = create_paper_level_data(df_with_cites)
    
    # Step 4: Calculate citation percentiles
    print("\\n" + "="*50)
    paper_with_flags = calculate_citation_percentiles(paper_data)
    
    # Step 5: Merge citation flags back to gene level
    print("\\n" + "="*50)
    df_enriched = merge_citation_flags_back(df_with_cites, paper_with_flags)
    
    # Step 6: Save results
    print("\\n" + "="*50)
    save_results(df_enriched, output_file, merge_stats)
    
    print(f"\\n=== Phase 4 Complete ===")
    print(f"Output: {output_file}")
    print(f"Citation-enriched records: {len(df_enriched):,}")
    print(f"With high-citation flags at annual/quarterly/bi-monthly levels")
    print(f"Ready for Phase 5: Panel Construction!")

if __name__ == "__main__":
    main()