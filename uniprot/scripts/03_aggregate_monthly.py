#!/usr/bin/env python3
"""
Aggregate publication data by month for temporal analysis.

This script takes the parsed publication data and creates monthly
time series datasets for analyzing research trends over time.

Usage:
    python 03_aggregate_monthly.py [--input-dir DATA_DIR] [--output-dir OUTPUT_DIR]

Author: AlphaFold Impact Analysis Pipeline
Date: October 2025
"""

import argparse
import logging
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('aggregate_monthly.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def create_monthly_time_series(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a monthly time series of publication counts per protein.

    Args:
        df: DataFrame with parsed publication data

    Returns:
        DataFrame with monthly counts per protein (long format)
    """
    logger.info("Creating monthly time series...")

    # Filter to only publications with month information
    df_with_month = df[df['month'].notna()].copy()
    logger.info(f"Using {len(df_with_month):,} publications with month information ({len(df_with_month)/len(df)*100:.1f}%)")

    # Create proper datetime for year-month
    df_with_month['year_month_date'] = pd.to_datetime(
        df_with_month[['year', 'month']].assign(day=1)
    )

    # Aggregate by protein and month
    monthly_counts = df_with_month.groupby(
        ['accession', 'gene_name', 'protein_name', 'year_month_date']
    ).agg({
        'pmid': 'count',  # Total publications
        'reviewed': lambda x: x.sum()  # Count of curated publications
    }).reset_index()

    monthly_counts.columns = ['accession', 'gene_name', 'protein_name', 'year_month', 'total_pubs', 'curated_pubs']

    # Calculate mapped publications
    monthly_counts['mapped_pubs'] = monthly_counts['total_pubs'] - monthly_counts['curated_pubs']

    # Sort by accession and date
    monthly_counts = monthly_counts.sort_values(['accession', 'year_month'])

    logger.info(f"Created time series with {len(monthly_counts):,} protein-month observations")

    return monthly_counts


def create_wide_format(monthly_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create wide format with months as columns for each protein.

    Args:
        monthly_df: Long format DataFrame with monthly counts

    Returns:
        Wide format DataFrame with one row per protein
    """
    logger.info("Creating wide format dataset...")

    # Pivot to wide format
    wide_df = monthly_df.pivot(
        index=['accession', 'gene_name', 'protein_name'],
        columns='year_month',
        values='total_pubs'
    ).reset_index()

    # Fill NaN with 0 (months with no publications)
    wide_df = wide_df.fillna(0)

    # Convert counts to integers
    month_columns = [col for col in wide_df.columns if isinstance(col, pd.Timestamp)]
    for col in month_columns:
        wide_df[col] = wide_df[col].astype(int)

    # Rename month columns to string format
    rename_dict = {col: col.strftime('%Y-%m') for col in month_columns}
    wide_df = wide_df.rename(columns=rename_dict)

    logger.info(f"Wide format: {len(wide_df):,} proteins × {len(month_columns)} months")

    return wide_df


def create_protein_summaries(df: pd.DataFrame, monthly_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary statistics for each protein.

    Args:
        df: Full parsed publication data
        monthly_df: Monthly aggregated data

    Returns:
        DataFrame with protein-level summaries
    """
    logger.info("Creating protein summaries...")

    # Get publications with month info
    df_with_month = df[df['month'].notna()].copy()
    df_with_month['year_month_date'] = pd.to_datetime(
        df_with_month[['year', 'month']].assign(day=1)
    )

    # Protein-level aggregations
    protein_stats = df.groupby(['accession', 'gene_name', 'protein_name']).agg({
        'pmid': 'count',  # Total publications
        'year': ['min', 'max'],  # First and last publication years
        'reviewed': 'first'  # Is this a reviewed protein?
    }).reset_index()

    protein_stats.columns = ['accession', 'gene_name', 'protein_name', 'total_pubs', 'first_pub_year', 'last_pub_year', 'is_reviewed']

    # Calculate stats from monthly data (only for pubs with month info)
    monthly_stats = df_with_month.groupby('accession').agg({
        'year_month_date': ['min', 'max'],
        'pmid': 'count'
    }).reset_index()
    monthly_stats.columns = ['accession', 'first_pub_date', 'last_pub_date', 'pubs_with_month']

    # Merge
    protein_stats = protein_stats.merge(monthly_stats, on='accession', how='left')

    # Calculate average monthly publications (for proteins with monthly data)
    monthly_avg = monthly_df.groupby('accession')['total_pubs'].mean().reset_index()
    monthly_avg.columns = ['accession', 'avg_monthly_pubs']
    protein_stats = protein_stats.merge(monthly_avg, on='accession', how='left')

    # Find peak publication month
    peak_months = monthly_df.loc[monthly_df.groupby('accession')['total_pubs'].idxmax()][['accession', 'year_month', 'total_pubs']]
    peak_months.columns = ['accession', 'peak_month', 'peak_month_pubs']
    protein_stats = protein_stats.merge(peak_months, on='accession', how='left')

    # Sort by total publications
    protein_stats = protein_stats.sort_values('total_pubs', ascending=False)

    logger.info(f"Created summaries for {len(protein_stats):,} proteins")

    return protein_stats


def create_yearly_aggregates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create yearly aggregates for broader temporal trends.

    Args:
        df: Full parsed publication data

    Returns:
        DataFrame with yearly counts per protein
    """
    logger.info("Creating yearly aggregates...")

    yearly_counts = df.groupby(['accession', 'gene_name', 'protein_name', 'year']).agg({
        'pmid': 'count',
        'reviewed': lambda x: x.sum()
    }).reset_index()

    yearly_counts.columns = ['accession', 'gene_name', 'protein_name', 'year', 'total_pubs', 'curated_pubs']
    yearly_counts['mapped_pubs'] = yearly_counts['total_pubs'] - yearly_counts['curated_pubs']

    yearly_counts = yearly_counts.sort_values(['accession', 'year'])

    logger.info(f"Created yearly aggregates: {len(yearly_counts):,} protein-year observations")

    return yearly_counts


def create_global_time_series(monthly_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create global time series (all proteins combined) by month.

    Args:
        monthly_df: Monthly data per protein

    Returns:
        DataFrame with total publications per month across all proteins
    """
    logger.info("Creating global monthly time series...")

    global_monthly = monthly_df.groupby('year_month').agg({
        'total_pubs': 'sum',
        'curated_pubs': 'sum',
        'mapped_pubs': 'sum',
        'accession': 'nunique'  # Number of proteins with publications that month
    }).reset_index()

    global_monthly.columns = ['year_month', 'total_pubs', 'curated_pubs', 'mapped_pubs', 'num_proteins']

    global_monthly = global_monthly.sort_values('year_month')

    logger.info(f"Created global time series: {len(global_monthly):,} months")

    return global_monthly


def generate_aggregation_stats(monthly_df: pd.DataFrame, protein_summaries: pd.DataFrame) -> Dict:
    """
    Generate statistics about the aggregated data.

    Args:
        monthly_df: Monthly aggregated data
        protein_summaries: Protein summary statistics

    Returns:
        Dictionary with aggregation statistics
    """
    stats = {
        'total_protein_months': len(monthly_df),
        'proteins_with_monthly_data': monthly_df['accession'].nunique(),
        'date_range': {
            'start': monthly_df['year_month'].min().strftime('%Y-%m'),
            'end': monthly_df['year_month'].max().strftime('%Y-%m')
        },
        'proteins_by_pub_count': {
            '0': int((protein_summaries['total_pubs'] == 0).sum()),
            '1-10': int(((protein_summaries['total_pubs'] >= 1) & (protein_summaries['total_pubs'] <= 10)).sum()),
            '11-50': int(((protein_summaries['total_pubs'] >= 11) & (protein_summaries['total_pubs'] <= 50)).sum()),
            '51-100': int(((protein_summaries['total_pubs'] >= 51) & (protein_summaries['total_pubs'] <= 100)).sum()),
            '101-500': int(((protein_summaries['total_pubs'] >= 101) & (protein_summaries['total_pubs'] <= 500)).sum()),
            '500+': int((protein_summaries['total_pubs'] > 500).sum())
        },
        'avg_monthly_pubs_per_protein': float(monthly_df.groupby('accession')['total_pubs'].mean().mean()),
        'median_monthly_pubs_per_protein': float(monthly_df.groupby('accession')['total_pubs'].mean().median())
    }

    return stats


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Aggregate publication data by month for temporal analysis"
    )
    parser.add_argument(
        '--input-dir',
        type=Path,
        default=Path('data/processed'),
        help='Input directory with parsed publication data (default: data/processed)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/outputs'),
        help='Output directory for aggregated datasets (default: data/outputs)'
    )

    args = parser.parse_args()

    logger.info("="*60)
    logger.info("UniProt Publication Monthly Aggregation")
    logger.info("="*60)
    logger.info(f"Input directory: {args.input_dir.absolute()}")
    logger.info(f"Output directory: {args.output_dir.absolute()}")
    logger.info("")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load parsed publication data
    input_file = args.input_dir / "all_publications.parquet"
    if not input_file.exists():
        input_file = args.input_dir / "all_publications.csv"

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        logger.error("Please run 02_parse_publications.py first")
        return

    logger.info(f"Loading data from {input_file}...")
    if input_file.suffix == '.parquet':
        df = pd.read_parquet(input_file)
    else:
        df = pd.read_csv(input_file)

    logger.info(f"Loaded {len(df):,} publications for {df['accession'].nunique():,} proteins")

    # Create monthly time series (long format)
    monthly_df = create_monthly_time_series(df)
    monthly_file = args.output_dir / "monthly_publications_long.csv"
    monthly_df.to_csv(monthly_file, index=False)
    logger.info(f"Saved monthly time series (long format) to {monthly_file}")

    # Create wide format
    wide_df = create_wide_format(monthly_df)
    wide_file = args.output_dir / "monthly_publications_wide.csv"
    wide_df.to_csv(wide_file, index=False)
    logger.info(f"Saved monthly time series (wide format) to {wide_file}")

    # Create protein summaries
    protein_summaries = create_protein_summaries(df, monthly_df)
    summary_file = args.output_dir / "protein_summaries.csv"
    protein_summaries.to_csv(summary_file, index=False)
    logger.info(f"Saved protein summaries to {summary_file}")

    # Create yearly aggregates
    yearly_df = create_yearly_aggregates(df)
    yearly_file = args.output_dir / "yearly_publications.csv"
    yearly_df.to_csv(yearly_file, index=False)
    logger.info(f"Saved yearly aggregates to {yearly_file}")

    # Create global time series
    global_monthly = create_global_time_series(monthly_df)
    global_file = args.output_dir / "global_monthly_timeseries.csv"
    global_monthly.to_csv(global_file, index=False)
    logger.info(f"Saved global monthly time series to {global_file}")

    # Generate and save statistics
    stats = generate_aggregation_stats(monthly_df, protein_summaries)
    stats_file = args.output_dir / "aggregation_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Saved aggregation statistics to {stats_file}")

    # Print summary
    logger.info("="*60)
    logger.info("AGGREGATION SUMMARY")
    logger.info("="*60)
    logger.info(f"Protein-month observations: {stats['total_protein_months']:,}")
    logger.info(f"Proteins with monthly data: {stats['proteins_with_monthly_data']:,}")
    logger.info(f"Date range: {stats['date_range']['start']} to {stats['date_range']['end']}")
    logger.info("")
    logger.info("Proteins by publication count:")
    for range_label, count in stats['proteins_by_pub_count'].items():
        logger.info(f"  {range_label}: {count:,} proteins")
    logger.info("")
    logger.info(f"Average monthly pubs per protein: {stats['avg_monthly_pubs_per_protein']:.2f}")
    logger.info(f"Median monthly pubs per protein: {stats['median_monthly_pubs_per_protein']:.2f}")

    logger.info("="*60)
    logger.info("Output files created:")
    logger.info("="*60)
    logger.info(f"  1. {monthly_file.name} - Monthly counts (long format)")
    logger.info(f"  2. {wide_file.name} - Monthly counts (wide format)")
    logger.info(f"  3. {summary_file.name} - Protein-level summaries")
    logger.info(f"  4. {yearly_file.name} - Yearly aggregates")
    logger.info(f"  5. {global_file.name} - Global monthly time series")
    logger.info(f"  6. {stats_file.name} - Aggregation statistics")
    logger.info("="*60)
    logger.info("Aggregation complete!")
    logger.info("="*60)


if __name__ == "__main__":
    main()
