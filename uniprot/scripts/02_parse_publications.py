#!/usr/bin/env python3
"""
Parse publication data from UniProt protein entries.

This script processes the downloaded UniProt JSON files and extracts
all publication references with their dates, creating a structured
dataset for temporal analysis.

Usage:
    python 02_parse_publications.py [--input-dir DATA_DIR] [--output-dir OUTPUT_DIR]

Author: AlphaFold Impact Analysis Pipeline
Date: October 2025
"""

import argparse
import logging
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('parse_publications.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def extract_publications_from_entry(protein_entry: dict, min_year: int = None, max_year: int = None) -> List[Dict]:
    """
    Extract all publication data from a single UniProt protein entry.

    Args:
        protein_entry: Dictionary containing a UniProt protein entry
        min_year: Minimum year to include (inclusive), None for no filter
        max_year: Maximum year to include (inclusive), None for no filter

    Returns:
        List of publication dictionaries with extracted fields
    """
    try:
        accession = protein_entry.get('primaryAccession', 'Unknown')

        # Extract gene name (primary)
        gene_name = 'Unknown'
        genes = protein_entry.get('genes', [])
        if genes and len(genes) > 0:
            gene_name = genes[0].get('geneName', {}).get('value', 'Unknown')

        # Extract protein name
        protein_name = protein_entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
        if protein_name == 'Unknown':
            # Try alternative name
            alt_names = protein_entry.get('proteinDescription', {}).get('alternativeNames', [])
            if alt_names:
                protein_name = alt_names[0].get('fullName', {}).get('value', 'Unknown')

        # Check if reviewed (Swiss-Prot vs TrEMBL)
        reviewed = protein_entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)'

        publications = []

        # Iterate through all references
        for ref in protein_entry.get('references', []):
            citation = ref.get('citation', {})

            # Extract publication date
            pub_date = citation.get('publicationDate')
            if not pub_date:
                continue  # Skip entries without dates

            # Parse date components
            date_parts = pub_date.split('-')
            year = None
            month = None

            try:
                if len(date_parts) >= 1:
                    year = int(date_parts[0])
                if len(date_parts) >= 2:
                    month = int(date_parts[1])
            except ValueError:
                logger.warning(f"Could not parse date '{pub_date}' for {accession}")
                continue

            # Apply year filter if specified
            if min_year is not None and year < min_year:
                continue
            if max_year is not None and year > max_year:
                continue

            # Extract PubMed ID
            pmid = None
            cross_refs = citation.get('citationCrossReferences', [])
            for xref in cross_refs:
                if xref.get('database') == 'PubMed':
                    pmid = xref.get('id')
                    break

            # Extract DOI
            doi = None
            for xref in cross_refs:
                if xref.get('database') == 'DOI':
                    doi = xref.get('id')
                    break

            # Extract citation type (curated vs mapped)
            citation_type = 'curated' if reviewed else 'mapped'

            # Extract reference categories (what this paper annotates)
            categories = []
            for scope in ref.get('referencePositions', []):
                categories.append(scope)

            # Extract journal information
            journal = citation.get('journal')
            if not journal:
                # For books or other citation types
                citation_type_name = citation.get('type', 'Unknown')
                journal = f"[{citation_type_name}]"

            # Build publication record
            pub_record = {
                'accession': accession,
                'gene_name': gene_name,
                'protein_name': protein_name,
                'reviewed': reviewed,
                'pmid': pmid,
                'doi': doi,
                'title': citation.get('title', ''),
                'publication_date': pub_date,
                'year': year,
                'month': month,
                'journal': journal,
                'citation_type': citation_type,
                'categories': '; '.join(categories) if categories else '',
                'authors': citation.get('authors', [])
            }

            publications.append(pub_record)

        return publications

    except Exception as e:
        logger.error(f"Error processing entry {protein_entry.get('primaryAccession', 'Unknown')}: {e}")
        return []


def parse_uniprot_file(input_file: Path, dataset_type: str = 'reviewed', min_year: int = None, max_year: int = None) -> pd.DataFrame:
    """
    Parse a UniProt JSON file and extract all publications.

    Args:
        input_file: Path to UniProt JSON file
        dataset_type: Type of dataset ('reviewed' or 'unreviewed')
        min_year: Minimum year to include (inclusive), None for no filter
        max_year: Maximum year to include (inclusive), None for no filter

    Returns:
        DataFrame containing all extracted publications
    """
    logger.info(f"Parsing {dataset_type} proteins from {input_file}")

    # Load JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)

    entries = data.get('results', [])
    logger.info(f"Found {len(entries):,} protein entries")

    if min_year or max_year:
        year_filter_msg = f"Filtering publications: "
        if min_year:
            year_filter_msg += f"{min_year} <= year"
        if max_year:
            if min_year:
                year_filter_msg += f" <= {max_year}"
            else:
                year_filter_msg += f"year <= {max_year}"
        logger.info(year_filter_msg)

    # Extract publications from all entries
    all_publications = []

    for entry in tqdm(entries, desc=f"Extracting {dataset_type} publications"):
        pubs = extract_publications_from_entry(entry, min_year=min_year, max_year=max_year)
        all_publications.extend(pubs)

    # Convert to DataFrame
    df = pd.DataFrame(all_publications)

    logger.info(f"Extracted {len(df):,} publications from {len(entries):,} proteins")

    return df


def clean_and_validate_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and validate the extracted publication data.

    Args:
        df: DataFrame with raw publication data

    Returns:
        Cleaned DataFrame
    """
    logger.info("Cleaning and validating data...")

    initial_count = len(df)

    # Remove entries with missing critical fields
    df = df.dropna(subset=['accession', 'year'])

    # Remove invalid years (e.g., future years or too old)
    current_year = datetime.now().year
    df = df[(df['year'] >= 1900) & (df['year'] <= current_year)]

    # Validate months (1-12)
    df.loc[df['month'].notna(), 'month'] = df.loc[df['month'].notna(), 'month'].astype(int)
    df = df[(df['month'].isna()) | ((df['month'] >= 1) & (df['month'] <= 12))]

    # Create a properly formatted date string
    def format_date(row):
        if pd.notna(row['month']):
            return f"{int(row['year'])}-{int(row['month']):02d}"
        else:
            return f"{int(row['year'])}"

    df['year_month'] = df.apply(format_date, axis=1)

    # Sort by accession and date
    df = df.sort_values(['accession', 'year', 'month'])

    removed_count = initial_count - len(df)
    logger.info(f"Removed {removed_count:,} invalid entries ({removed_count/initial_count*100:.1f}%)")
    logger.info(f"Final dataset: {len(df):,} publications")

    return df


def generate_summary_statistics(df: pd.DataFrame) -> Dict:
    """
    Generate summary statistics about the parsed publication data.

    Args:
        df: DataFrame with publication data

    Returns:
        Dictionary with summary statistics
    """
    stats = {
        'total_publications': len(df),
        'total_proteins': df['accession'].nunique(),
        'total_genes': df['gene_name'].nunique(),
        'reviewed_proteins': df[df['reviewed']]['accession'].nunique(),
        'unreviewed_proteins': df[~df['reviewed']]['accession'].nunique(),
        'publications_with_pmid': df['pmid'].notna().sum(),
        'publications_with_doi': df['doi'].notna().sum(),
        'publications_with_month': df['month'].notna().sum(),
        'year_range': {
            'min': int(df['year'].min()),
            'max': int(df['year'].max())
        },
        'top_10_proteins_by_pub_count': df.groupby(['accession', 'gene_name']).size().sort_values(ascending=False).head(10).to_dict(),
        'publications_per_year': df.groupby('year').size().to_dict(),
        'avg_pubs_per_protein': float(df.groupby('accession').size().mean()),
        'median_pubs_per_protein': float(df.groupby('accession').size().median())
    }

    return stats


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Parse publication data from UniProt protein entries"
    )
    parser.add_argument(
        '--input-dir',
        type=Path,
        default=Path('data/raw'),
        help='Input directory with downloaded UniProt files (default: data/raw)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/processed'),
        help='Output directory for parsed data (default: data/processed)'
    )
    parser.add_argument(
        '--min-year',
        type=int,
        default=2020,
        help='Minimum publication year to include (default: 2020)'
    )
    parser.add_argument(
        '--max-year',
        type=int,
        default=2023,
        help='Maximum publication year to include (default: 2023)'
    )

    args = parser.parse_args()

    logger.info("="*60)
    logger.info("UniProt Publication Data Parser")
    logger.info("="*60)
    logger.info(f"Input directory: {args.input_dir.absolute()}")
    logger.info(f"Output directory: {args.output_dir.absolute()}")
    logger.info(f"Year range: {args.min_year} - {args.max_year}")
    logger.info("")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_publications = []

    # Parse reviewed proteins
    reviewed_file = args.input_dir / "human_reviewed.json"
    if reviewed_file.exists():
        logger.info("Processing reviewed (Swiss-Prot) proteins...")
        reviewed_df = parse_uniprot_file(reviewed_file, 'reviewed', min_year=args.min_year, max_year=args.max_year)
        all_publications.append(reviewed_df)
    else:
        logger.warning(f"Reviewed proteins file not found: {reviewed_file}")

    # Parse unreviewed proteins
    unreviewed_file = args.input_dir / "human_unreviewed.json"
    if unreviewed_file.exists():
        logger.info("Processing unreviewed (TrEMBL) proteins...")
        unreviewed_df = parse_uniprot_file(unreviewed_file, 'unreviewed', min_year=args.min_year, max_year=args.max_year)
        all_publications.append(unreviewed_df)
    else:
        logger.warning(f"Unreviewed proteins file not found: {unreviewed_file}")

    if not all_publications:
        logger.error("No data files found to process!")
        return

    # Combine all publications
    logger.info("Combining all publications...")
    combined_df = pd.concat(all_publications, ignore_index=True)

    # Clean and validate
    cleaned_df = clean_and_validate_data(combined_df)

    # Save full dataset
    output_file = args.output_dir / "all_publications.csv"
    cleaned_df.to_csv(output_file, index=False)
    logger.info(f"Saved full publication dataset to {output_file}")

    # Save in Parquet format for efficient storage
    parquet_file = args.output_dir / "all_publications.parquet"
    cleaned_df.to_parquet(parquet_file, index=False, compression='snappy')
    logger.info(f"Saved Parquet format to {parquet_file}")

    # Generate and save summary statistics
    stats = generate_summary_statistics(cleaned_df)
    stats_file = args.output_dir / "parsing_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Saved summary statistics to {stats_file}")

    # Print summary
    logger.info("="*60)
    logger.info("PARSING SUMMARY")
    logger.info("="*60)
    logger.info(f"Total publications: {stats['total_publications']:,}")
    logger.info(f"Total proteins: {stats['total_proteins']:,}")
    logger.info(f"  - Reviewed: {stats['reviewed_proteins']:,}")
    logger.info(f"  - Unreviewed: {stats['unreviewed_proteins']:,}")
    logger.info(f"Publications with PMID: {stats['publications_with_pmid']:,} ({stats['publications_with_pmid']/stats['total_publications']*100:.1f}%)")
    logger.info(f"Publications with month: {stats['publications_with_month']:,} ({stats['publications_with_month']/stats['total_publications']*100:.1f}%)")
    logger.info(f"Year range: {stats['year_range']['min']} - {stats['year_range']['max']}")
    logger.info(f"Average pubs per protein: {stats['avg_pubs_per_protein']:.1f}")
    logger.info(f"Median pubs per protein: {stats['median_pubs_per_protein']:.1f}")

    logger.info("\nTop 10 proteins by publication count:")
    for (accession, gene), count in list(stats['top_10_proteins_by_pub_count'].items())[:10]:
        logger.info(f"  {gene} ({accession}): {count:,} publications")

    logger.info("="*60)
    logger.info("Parsing complete!")
    logger.info("="*60)


if __name__ == "__main__":
    main()
