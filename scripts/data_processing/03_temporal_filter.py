#!/usr/bin/env python3
"""
Phase 3: Temporal Filtering

Add publication dates and filter to 2020-2023 time window for AlphaFold impact analysis.

Input:  
- processed/gene_disease_intersect.parquet (33.6M rows, 118 MB)
- raw_data/complete_pmid_dates.csv (6.7M PMIDs with dates)

Output: processed/gene_disease_temporal.parquet

Processing steps:
1. Load gene-disease intersection data from Phase 2
2. Merge with publication dates from NIH data  
3. Filter to 2020-2023 time window (AlphaFold era)
4. Create time variables (year, month, ym, quarters, bi-monthly bins)
5. Document temporal distribution and merge success
6. Save with full statistics
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_gene_disease_data(file_path):
    """Load gene-disease intersection data from Phase 2."""
    print(f"Loading gene-disease intersection data from: {file_path}")
    
    df = pd.read_parquet(file_path)
    
    print(f"  Loaded {len(df):,} gene-disease intersection records")
    print(f"  Unique PMIDs: {df['pmid'].nunique():,}")
    print(f"  Unique genes: {df['gene_id'].nunique():,}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def load_publication_dates(file_path):
    """Load publication dates from NIH data."""
    print(f"Loading publication dates from: {file_path}")
    
    # Read without dtype constraints first to handle missing values
    df = pd.read_csv(file_path)
    
    initial_count = len(df)
    print(f"  Loaded {initial_count:,} PMID-date records")
    
    # Data quality checks before cleaning
    null_years = df['year'].isnull().sum()
    null_months = df['month'].isnull().sum()
    
    if null_years > 0 or null_months > 0:
        print(f"  Found missing values - Year: {null_years:,}, Month: {null_months:,}")
        
        # Drop records with missing dates
        df = df.dropna(subset=['year', 'month']).copy()
        print(f"  After dropping missing dates: {len(df):,} records ({len(df)/initial_count*100:.2f}% retained)")
    
    # Convert to appropriate dtypes
    df['pmid'] = df['pmid'].astype('int32')
    df['year'] = df['year'].astype('int16') 
    df['month'] = df['month'].astype('int8')
    
    print(f"  Date range: {df['year'].min()}-{df['year'].max()}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    # Additional data quality checks
    invalid_months = (~df['month'].between(1, 12)).sum()
    if invalid_months > 0:
        print(f"  ⚠️ Warning: {invalid_months:,} records with invalid months")
        df = df[df['month'].between(1, 12)].copy()
        print(f"  After filtering invalid months: {len(df):,} records")
    
    return df

def merge_dates(df_gene, df_dates):
    """Merge gene-disease data with publication dates."""
    print("Merging gene-disease data with publication dates...")
    
    initial_count = len(df_gene)
    initial_pmids = df_gene['pmid'].nunique()
    
    # Merge on pmid
    df_merged = df_gene.merge(df_dates, on='pmid', how='inner')
    
    final_count = len(df_merged)
    final_pmids = df_merged['pmid'].nunique()
    
    merge_rate = final_count / initial_count
    pmid_merge_rate = final_pmids / initial_pmids
    
    print(f"  Initial records: {initial_count:,}")
    print(f"  Merged records: {final_count:,} ({merge_rate*100:.2f}% retained)")
    print(f"  Initial PMIDs: {initial_pmids:,}")
    print(f"  Merged PMIDs: {final_pmids:,} ({pmid_merge_rate*100:.2f}% retained)")
    
    if merge_rate < 0.9:
        print(f"  ⚠️ Warning: Low merge rate - many PMIDs lack date information")
        
        # Sample of unmatched PMIDs for investigation
        unmatched_pmids = df_gene[~df_gene['pmid'].isin(df_dates['pmid'])]['pmid'].head(5)
        print(f"  Sample unmatched PMIDs: {unmatched_pmids.tolist()}")
    
    return df_merged, {
        'initial_records': initial_count,
        'final_records': final_count,
        'initial_pmids': initial_pmids,  
        'final_pmids': final_pmids,
        'merge_rate': merge_rate,
        'pmid_merge_rate': pmid_merge_rate
    }

def create_time_variables(df):
    """Create comprehensive time variables for analysis."""
    print("Creating time variables...")
    
    # Basic time variables
    df = df.copy()
    
    # Year-month combined (YYYYMM format for easy sorting)
    df['ym'] = df['year'] * 100 + df['month']
    
    # Quarterly variables
    df['quarter'] = ((df['month'] - 1) // 3) + 1  # 1-4 quarters per year
    df['year_quarter'] = df['year'] * 10 + df['quarter']  # YYYYQ format
    
    # Bi-monthly variables (6 periods per year)
    df['bimonth'] = ((df['month'] - 1) // 2) + 1  # 1-6 bi-months per year
    df['year_bimonth'] = df['year'] * 10 + df['bimonth']  # YYYYB format
    
    print(f"  Created time variables:")
    print(f"  Year range: {df['year'].min()}-{df['year'].max()}")
    print(f"  Month range: {df['month'].min()}-{df['month'].max()}")
    print(f"  YM range: {df['ym'].min()}-{df['ym'].max()}")
    
    return df

def filter_alphafold_era(df, start_year=2020, end_year=2023):
    """Filter to AlphaFold era (2020-2023) for impact analysis."""
    print(f"Filtering to AlphaFold era: {start_year}-{end_year}...")
    
    initial_count = len(df)
    initial_pmids = df['pmid'].nunique()
    initial_genes = df['gene_id'].nunique()
    
    # Show temporal distribution before filtering
    year_dist = df['year'].value_counts().sort_index()
    print(f"  Pre-filter year distribution:")
    for year, count in year_dist.items():
        print(f"    {year}: {count:,} records ({count/initial_count*100:.1f}%)")
    
    # Filter to target years
    df_filtered = df[df['year'].between(start_year, end_year)].copy()
    
    final_count = len(df_filtered)
    final_pmids = df_filtered['pmid'].nunique()
    final_genes = df_filtered['gene_id'].nunique()
    
    print(f"  Initial records: {initial_count:,}")
    print(f"  Filtered records: {final_count:,} ({final_count/initial_count*100:.2f}% retained)")
    print(f"  Initial PMIDs: {initial_pmids:,}")
    print(f"  Filtered PMIDs: {final_pmids:,} ({final_pmids/initial_pmids*100:.2f}% retained)")
    print(f"  Gene retention: {final_genes:,}/{initial_genes:,} ({final_genes/initial_genes*100:.2f}%)")
    
    # Post-filter year distribution
    year_dist_final = df_filtered['year'].value_counts().sort_index()
    print(f"  Post-filter year distribution:")
    for year, count in year_dist_final.items():
        print(f"    {year}: {count:,} records ({count/final_count*100:.1f}%)")
    
    return df_filtered, {
        'initial_records': initial_count,
        'final_records': final_count,
        'initial_pmids': initial_pmids,
        'final_pmids': final_pmids,
        'initial_genes': initial_genes,
        'final_genes': final_genes,
        'retention_rate': final_count/initial_count,
        'year_distribution': year_dist_final.to_dict()
    }

def save_results(df, output_path, merge_stats, filter_stats):
    """Save results with comprehensive metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} temporal records")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase3_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 3: Temporal Filtering Processing Statistics\n")
        f.write("=" * 55 + "\n\n")
        
        f.write("Input files:\n")
        f.write("  - processed/gene_disease_intersect.parquet (Phase 2 output)\n")
        f.write("  - raw_data/complete_pmid_dates.csv (NIH publication dates)\n\n")
        
        f.write(f"Output file: {output_path.name}\n\n")
        
        f.write("STEP 1: Date Merging\n")
        f.write("-" * 20 + "\n")
        f.write(f"Initial records: {merge_stats['initial_records']:,}\n")
        f.write(f"Records with dates: {merge_stats['final_records']:,}\n")
        f.write(f"Merge success rate: {merge_stats['merge_rate']*100:.2f}%\n")
        f.write(f"PMID merge rate: {merge_stats['pmid_merge_rate']*100:.2f}%\n\n")
        
        f.write("STEP 2: AlphaFold Era Filtering (2020-2023)\n")
        f.write("-" * 40 + "\n")
        f.write(f"Pre-filter records: {filter_stats['initial_records']:,}\n")
        f.write(f"Post-filter records: {filter_stats['final_records']:,}\n")
        f.write(f"Temporal retention: {filter_stats['retention_rate']*100:.2f}%\n\n")
        
        f.write(f"Pre-filter PMIDs: {filter_stats['initial_pmids']:,}\n")
        f.write(f"Post-filter PMIDs: {filter_stats['final_pmids']:,}\n\n")
        
        f.write(f"Pre-filter genes: {filter_stats['initial_genes']:,}\n")
        f.write(f"Post-filter genes: {filter_stats['final_genes']:,}\n\n")
        
        f.write("Final Year Distribution:\n")
        for year, count in sorted(filter_stats['year_distribution'].items()):
            pct = count / filter_stats['final_records'] * 100
            f.write(f"  {year}: {count:,} records ({pct:.1f}%)\n")
        f.write("\n")
        
        f.write("Final Output Summary\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Unique PMIDs: {df['pmid'].nunique():,}\n")
        f.write(f"Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"Time span: 2020-2023 (AlphaFold era)\n")
        f.write(f"Avg mentions per paper: {len(df)/df['pmid'].nunique():.2f}\n")
        f.write(f"Output file size: {file_size:.1f} MB\n\n")
        
        f.write("Time Variables Created:\n")
        f.write("  - year, month: Basic temporal units\n")
        f.write("  - ym: Year-month combined (YYYYMM)\n")
        f.write("  - quarter, year_quarter: Quarterly aggregation\n")
        f.write("  - bimonth, year_bimonth: Bi-monthly aggregation\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    gene_disease_input = base_dir / "processed" / "gene_disease_intersect.parquet"
    dates_input = base_dir / "raw_data" / "complete_pmid_dates.csv"
    output_file = base_dir / "processed" / "gene_disease_temporal.parquet"
    
    print("=== Phase 3: Temporal Filtering ===\n")
    
    # Step 1: Load gene-disease intersection data
    df_gene = load_gene_disease_data(gene_disease_input)
    
    # Step 2: Load publication dates  
    df_dates = load_publication_dates(dates_input)
    
    # Step 3: Merge with dates
    df_with_dates, merge_stats = merge_dates(df_gene, df_dates)
    
    # Step 4: Create time variables
    df_with_time = create_time_variables(df_with_dates)
    
    # Step 5: Filter to AlphaFold era (2020-2023)
    df_filtered, filter_stats = filter_alphafold_era(df_with_time)
    
    # Step 6: Save results
    save_results(df_filtered, output_file, merge_stats, filter_stats)
    
    print(f"\n=== Phase 3 Complete ===")
    print(f"Output: {output_file}")
    print(f"AlphaFold era records: {len(df_filtered):,}")
    print(f"Time span: 2020-2023")
    print(f"Unique genes: {df_filtered['gene_id'].nunique():,}")
    print(f"Unique PMIDs: {df_filtered['pmid'].nunique():,}")

if __name__ == "__main__":
    main()