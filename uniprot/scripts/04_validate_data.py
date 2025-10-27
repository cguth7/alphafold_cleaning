#!/usr/bin/env python3
"""
Validate and perform quality checks on extracted publication data.

This script performs comprehensive validation of the parsed and aggregated
publication data to ensure data quality and identify potential issues.

Usage:
    python 04_validate_data.py [--data-dir DATA_DIR]

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
from typing import Dict, List, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('validate_data.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class DataValidator:
    """Validator for UniProt publication data quality checks."""

    def __init__(self, processed_dir: Path, outputs_dir: Path):
        """
        Initialize validator with data directories.

        Args:
            processed_dir: Directory with processed data
            outputs_dir: Directory with aggregated outputs
        """
        self.processed_dir = processed_dir
        self.outputs_dir = outputs_dir
        self.validation_results = {
            'timestamp': datetime.now().isoformat(),
            'checks_passed': 0,
            'checks_failed': 0,
            'warnings': 0,
            'errors': []
        }

    def check_file_exists(self, file_path: Path, description: str) -> bool:
        """Check if a required file exists."""
        if file_path.exists():
            logger.info(f"✓ {description}: Found at {file_path}")
            self.validation_results['checks_passed'] += 1
            return True
        else:
            logger.error(f"✗ {description}: NOT FOUND at {file_path}")
            self.validation_results['checks_failed'] += 1
            self.validation_results['errors'].append(f"Missing file: {file_path}")
            return False

    def validate_publications_data(self, df: pd.DataFrame) -> Dict:
        """
        Validate the parsed publications dataset.

        Args:
            df: DataFrame with parsed publication data

        Returns:
            Dictionary with validation results
        """
        logger.info("="*60)
        logger.info("Validating Publications Data")
        logger.info("="*60)

        results = {}

        # Check 1: Required columns
        required_cols = ['accession', 'gene_name', 'year', 'publication_date']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"✗ Missing required columns: {missing_cols}")
            self.validation_results['checks_failed'] += 1
            results['missing_columns'] = missing_cols
        else:
            logger.info(f"✓ All required columns present")
            self.validation_results['checks_passed'] += 1

        # Check 2: Data completeness
        logger.info("\nData Completeness:")
        for col in required_cols:
            if col in df.columns:
                null_count = df[col].isna().sum()
                null_pct = (null_count / len(df)) * 100
                if null_pct > 5:  # Warn if >5% missing
                    logger.warning(f"  ⚠ {col}: {null_count:,} null values ({null_pct:.1f}%)")
                    self.validation_results['warnings'] += 1
                else:
                    logger.info(f"  ✓ {col}: {null_count:,} null values ({null_pct:.1f}%)")
                    self.validation_results['checks_passed'] += 1

        results['null_counts'] = df.isnull().sum().to_dict()

        # Check 3: Year validity
        current_year = datetime.now().year
        invalid_years = df[(df['year'] < 1900) | (df['year'] > current_year)]
        if len(invalid_years) > 0:
            logger.error(f"✗ Found {len(invalid_years):,} publications with invalid years")
            self.validation_results['checks_failed'] += 1
            results['invalid_years'] = len(invalid_years)
        else:
            logger.info(f"✓ All years valid (1900-{current_year})")
            self.validation_results['checks_passed'] += 1

        # Check 4: Month validity
        if 'month' in df.columns:
            invalid_months = df[df['month'].notna() & ((df['month'] < 1) | (df['month'] > 12))]
            if len(invalid_months) > 0:
                logger.error(f"✗ Found {len(invalid_months):,} publications with invalid months")
                self.validation_results['checks_failed'] += 1
                results['invalid_months'] = len(invalid_months)
            else:
                logger.info(f"✓ All months valid (1-12)")
                self.validation_results['checks_passed'] += 1

        # Check 5: Accession format
        invalid_accessions = df[~df['accession'].str.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', na=False)]
        if len(invalid_accessions) > 0:
            logger.warning(f"⚠ Found {len(invalid_accessions):,} entries with non-standard accession format")
            self.validation_results['warnings'] += 1
            results['non_standard_accessions'] = len(invalid_accessions)
        else:
            logger.info(f"✓ All accessions follow UniProt format")
            self.validation_results['checks_passed'] += 1

        # Check 6: Data distribution checks
        logger.info("\nData Distribution:")
        proteins_count = df['accession'].nunique()
        logger.info(f"  Total publications: {len(df):,}")
        logger.info(f"  Unique proteins: {proteins_count:,}")
        logger.info(f"  Avg pubs per protein: {len(df)/proteins_count:.1f}")

        pubs_per_protein = df.groupby('accession').size()
        results['distribution'] = {
            'mean': float(pubs_per_protein.mean()),
            'median': float(pubs_per_protein.median()),
            'std': float(pubs_per_protein.std()),
            'min': int(pubs_per_protein.min()),
            'max': int(pubs_per_protein.max())
        }

        # Check 7: Well-studied proteins
        logger.info("\nTop 5 most-studied proteins (validation check):")
        top_proteins = pubs_per_protein.sort_values(ascending=False).head(5)
        for acc, count in top_proteins.items():
            gene_name = df[df['accession'] == acc]['gene_name'].iloc[0]
            logger.info(f"  {gene_name} ({acc}): {count:,} publications")

        # Sanity check: Top proteins should have >100 publications for human proteins
        if top_proteins.iloc[0] < 100:
            logger.warning(f"⚠ Top protein has only {top_proteins.iloc[0]} publications - seems low for human proteome")
            self.validation_results['warnings'] += 1
        else:
            logger.info(f"✓ Top protein has {top_proteins.iloc[0]} publications - reasonable")
            self.validation_results['checks_passed'] += 1

        return results

    def validate_monthly_data(self, monthly_df: pd.DataFrame) -> Dict:
        """
        Validate monthly aggregated data.

        Args:
            monthly_df: DataFrame with monthly publication counts

        Returns:
            Dictionary with validation results
        """
        logger.info("="*60)
        logger.info("Validating Monthly Aggregated Data")
        logger.info("="*60)

        results = {}

        # Check 1: Date continuity
        if 'year_month' in monthly_df.columns:
            monthly_df['year_month'] = pd.to_datetime(monthly_df['year_month'])
            date_range = pd.date_range(
                start=monthly_df['year_month'].min(),
                end=monthly_df['year_month'].max(),
                freq='MS'
            )
            expected_months = len(date_range)
            actual_months = monthly_df['year_month'].nunique()

            logger.info(f"  Date range: {monthly_df['year_month'].min().strftime('%Y-%m')} to {monthly_df['year_month'].max().strftime('%Y-%m')}")
            logger.info(f"  Expected months: {expected_months}")
            logger.info(f"  Months with data: {actual_months}")

            # Not all months will have data for all proteins, so this is informational
            results['date_coverage'] = {
                'start': monthly_df['year_month'].min().isoformat(),
                'end': monthly_df['year_month'].max().isoformat(),
                'expected_months': expected_months,
                'months_with_data': actual_months
            }

        # Check 2: Negative counts
        count_cols = [col for col in monthly_df.columns if 'pubs' in col or 'count' in col]
        for col in count_cols:
            if col in monthly_df.columns and pd.api.types.is_numeric_dtype(monthly_df[col]):
                negative_count = (monthly_df[col] < 0).sum()
                if negative_count > 0:
                    logger.error(f"✗ {col}: Found {negative_count} negative values!")
                    self.validation_results['checks_failed'] += 1
                else:
                    logger.info(f"✓ {col}: No negative values")
                    self.validation_results['checks_passed'] += 1

        # Check 3: Curated vs mapped consistency
        if 'total_pubs' in monthly_df.columns and 'curated_pubs' in monthly_df.columns and 'mapped_pubs' in monthly_df.columns:
            inconsistent = monthly_df[monthly_df['total_pubs'] != monthly_df['curated_pubs'] + monthly_df['mapped_pubs']]
            if len(inconsistent) > 0:
                logger.error(f"✗ Found {len(inconsistent)} rows where total ≠ curated + mapped")
                self.validation_results['checks_failed'] += 1
            else:
                logger.info(f"✓ Curated + mapped = total for all rows")
                self.validation_results['checks_passed'] += 1

        return results

    def validate_protein_summaries(self, summary_df: pd.DataFrame) -> Dict:
        """
        Validate protein summary statistics.

        Args:
            summary_df: DataFrame with protein summaries

        Returns:
            Dictionary with validation results
        """
        logger.info("="*60)
        logger.info("Validating Protein Summaries")
        logger.info("="*60)

        results = {}

        # Check 1: Date logic
        if 'first_pub_year' in summary_df.columns and 'last_pub_year' in summary_df.columns:
            invalid_dates = summary_df[summary_df['first_pub_year'] > summary_df['last_pub_year']]
            if len(invalid_dates) > 0:
                logger.error(f"✗ Found {len(invalid_dates)} proteins where first_pub_year > last_pub_year")
                self.validation_results['checks_failed'] += 1
            else:
                logger.info(f"✓ All first publication dates ≤ last publication dates")
                self.validation_results['checks_passed'] += 1

        # Check 2: Publication counts
        if 'total_pubs' in summary_df.columns:
            zero_pubs = (summary_df['total_pubs'] == 0).sum()
            if zero_pubs > 0:
                logger.info(f"  ℹ {zero_pubs:,} proteins with 0 publications (expected for uncharacterized proteins)")

            # Distribution
            logger.info(f"\nPublication count distribution:")
            logger.info(f"  Mean: {summary_df['total_pubs'].mean():.1f}")
            logger.info(f"  Median: {summary_df['total_pubs'].median():.1f}")
            logger.info(f"  Max: {summary_df['total_pubs'].max():,}")

            results['pub_count_distribution'] = {
                'mean': float(summary_df['total_pubs'].mean()),
                'median': float(summary_df['total_pubs'].median()),
                'max': int(summary_df['total_pubs'].max()),
                'proteins_with_zero': int(zero_pubs)
            }

        # Check 3: Expected major proteins
        expected_proteins = ['TP53', 'BRCA1', 'EGFR', 'TNF', 'INS', 'ALB']
        found_proteins = summary_df[summary_df['gene_name'].isin(expected_proteins)]

        logger.info(f"\nExpected well-studied proteins:")
        for _, row in found_proteins.iterrows():
            logger.info(f"  ✓ {row['gene_name']} ({row['accession']}): {row['total_pubs']:,} publications")

        missing = set(expected_proteins) - set(found_proteins['gene_name'])
        if missing:
            logger.warning(f"  ⚠ Missing expected proteins: {missing}")
            self.validation_results['warnings'] += 1

        return results

    def validate_temporal_consistency(self, monthly_df: pd.DataFrame, yearly_df: pd.DataFrame) -> Dict:
        """
        Validate consistency between monthly and yearly aggregates.

        Args:
            monthly_df: Monthly aggregated data
            yearly_df: Yearly aggregated data

        Returns:
            Dictionary with validation results
        """
        logger.info("="*60)
        logger.info("Validating Temporal Consistency")
        logger.info("="*60)

        results = {}

        # Sum monthly data to yearly level
        if 'year_month' in monthly_df.columns:
            monthly_df['year'] = pd.to_datetime(monthly_df['year_month']).dt.year
            monthly_to_yearly = monthly_df.groupby(['accession', 'year'])['total_pubs'].sum().reset_index()

            # Compare with yearly data
            merged = monthly_to_yearly.merge(
                yearly_df[['accession', 'year', 'total_pubs']],
                on=['accession', 'year'],
                suffixes=('_monthly', '_yearly'),
                how='outer'
            )

            mismatches = merged[merged['total_pubs_monthly'] != merged['total_pubs_yearly']].dropna()

            if len(mismatches) > 0:
                logger.warning(f"⚠ Found {len(mismatches)} year-protein combinations with monthly/yearly mismatches")
                logger.warning(f"  This is expected if some publications lack month information")
                self.validation_results['warnings'] += 1
                results['monthly_yearly_mismatches'] = len(mismatches)
            else:
                logger.info(f"✓ Monthly aggregations match yearly data")
                self.validation_results['checks_passed'] += 1

        return results

    def run_all_validations(self) -> Dict:
        """
        Run all validation checks.

        Returns:
            Complete validation results dictionary
        """
        logger.info("="*60)
        logger.info("UniProt Data Validation Suite")
        logger.info("="*60)

        # Check file existence
        pub_file = self.processed_dir / "all_publications.parquet"
        if not pub_file.exists():
            pub_file = self.processed_dir / "all_publications.csv"

        self.check_file_exists(pub_file, "Parsed publications dataset")
        self.check_file_exists(self.outputs_dir / "monthly_publications_long.csv", "Monthly publications (long format)")
        self.check_file_exists(self.outputs_dir / "monthly_publications_wide.csv", "Monthly publications (wide format)")
        self.check_file_exists(self.outputs_dir / "protein_summaries.csv", "Protein summaries")
        self.check_file_exists(self.outputs_dir / "yearly_publications.csv", "Yearly publications")
        self.check_file_exists(self.outputs_dir / "global_monthly_timeseries.csv", "Global monthly time series")

        # Load data
        logger.info("\nLoading datasets for validation...")
        try:
            if pub_file.suffix == '.parquet':
                publications_df = pd.read_parquet(pub_file)
            else:
                publications_df = pd.read_csv(pub_file)

            monthly_df = pd.read_csv(self.outputs_dir / "monthly_publications_long.csv")
            yearly_df = pd.read_csv(self.outputs_dir / "yearly_publications.csv")
            summary_df = pd.read_csv(self.outputs_dir / "protein_summaries.csv")

            # Run validation checks
            self.validation_results['publications'] = self.validate_publications_data(publications_df)
            self.validation_results['monthly'] = self.validate_monthly_data(monthly_df)
            self.validation_results['summaries'] = self.validate_protein_summaries(summary_df)
            self.validation_results['temporal_consistency'] = self.validate_temporal_consistency(monthly_df, yearly_df)

        except Exception as e:
            logger.error(f"Error during validation: {e}")
            self.validation_results['errors'].append(str(e))
            self.validation_results['checks_failed'] += 1

        return self.validation_results


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Validate UniProt publication data quality"
    )
    parser.add_argument(
        '--processed-dir',
        type=Path,
        default=Path('data/processed'),
        help='Directory with processed data (default: data/processed)'
    )
    parser.add_argument(
        '--outputs-dir',
        type=Path,
        default=Path('data/outputs'),
        help='Directory with aggregated outputs (default: data/outputs)'
    )
    parser.add_argument(
        '--output-file',
        type=Path,
        default=Path('data/outputs/validation_results.json'),
        help='Output file for validation results (default: data/outputs/validation_results.json)'
    )

    args = parser.parse_args()

    # Create validator
    validator = DataValidator(args.processed_dir, args.outputs_dir)

    # Run validations
    results = validator.run_all_validations()

    # Save results
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_file, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    logger.info("="*60)
    logger.info("VALIDATION SUMMARY")
    logger.info("="*60)
    logger.info(f"Checks passed: {results['checks_passed']}")
    logger.info(f"Checks failed: {results['checks_failed']}")
    logger.info(f"Warnings: {results['warnings']}")

    if results['errors']:
        logger.info(f"\nErrors encountered:")
        for error in results['errors']:
            logger.error(f"  - {error}")

    logger.info(f"\nFull validation results saved to: {args.output_file}")

    # Exit code based on failures
    if results['checks_failed'] > 0:
        logger.error("\n⚠ VALIDATION FAILED - Please review errors above")
        return 1
    else:
        logger.info("\n✓ ALL VALIDATIONS PASSED")
        return 0


if __name__ == "__main__":
    exit(main())
