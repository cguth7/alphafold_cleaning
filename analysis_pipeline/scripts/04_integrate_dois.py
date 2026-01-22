"""
04_integrate_dois.py - Integrate DOI data into the analysis

This script:
1. Loads the DOI-PMID crosswalk
2. Creates a gene-semester level DOI mapping for downstream analysis
3. Exports DOI list for Danilo's Shi & Evans metrics

Output:
    - data/gene_semester_dois.parquet: Gene-semester aggregated DOI data
    - data/dois_for_shi_evans.csv: DOI list for external analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "analysis_pipeline" / "data"
DOI_FILE = BASE_DIR / "new_stuff_for_claude" / "doi_pmid_date_gene_2017_2025.csv"

# Treatment timing
TREATMENT_YM = 738  # July 2021


def load_doi_crosswalk():
    """Load the DOI-PMID-gene crosswalk."""
    print(f"Loading DOI crosswalk from {DOI_FILE}")

    # Load in chunks to handle large file
    chunks = []
    for chunk in pd.read_csv(DOI_FILE, chunksize=1_000_000):
        chunks.append(chunk)

    df = pd.concat(chunks, ignore_index=True)
    print(f"  Loaded {len(df):,} records")
    print(f"  Columns: {df.columns.tolist()}")
    print(f"  Year range: {df['year'].min()} - {df['year'].max()}")

    return df


def create_time_variables(df):
    """Add time variables matching the main analysis."""
    df = df.copy()

    # Drop rows with missing year/month
    df = df.dropna(subset=['year', 'month'])

    # Create ym (Stata-style year-month)
    df['ym'] = (df['year'] - 1960) * 12 + (df['month'] - 1)

    # Semester
    ym_min = 720  # Jan 2020
    df['ym_seq'] = df['ym'] - ym_min + 1
    df['semester'] = ((df['ym_seq'] - 1) // 6).astype(int)

    # Relative semester
    treatment_seq = TREATMENT_YM - ym_min + 1
    treatment_semester = (treatment_seq - 1) // 6
    df['rel_semester'] = df['semester'] - treatment_semester

    return df


def aggregate_dois_by_gene_semester(df):
    """Aggregate DOIs to gene-semester level."""
    print("Aggregating DOIs by gene-semester...")

    # Filter to analysis period (2020-2023)
    df = df[(df['year'] >= 2020) & (df['year'] <= 2023)].copy()

    # Add time variables
    df = create_time_variables(df)

    # Aggregate: list of DOIs per gene-semester
    agg = df.groupby(['gene_id', 'semester', 'rel_semester']).agg({
        'doi': list,
        'pmid': 'count'
    }).reset_index()

    agg.columns = ['gene_id', 'semester', 'rel_semester', 'dois', 'n_papers_doi']

    print(f"  Created {len(agg):,} gene-semester records")

    return agg


def export_dois_for_shi_evans(df):
    """Export DOI list for Shi & Evans analysis."""
    print("Exporting DOIs for Shi & Evans...")

    # Get unique DOIs with metadata
    dois = df[['doi', 'pmid', 'year', 'month', 'gene_id']].drop_duplicates(subset=['doi'])

    # Filter to analysis period
    dois = dois[(dois['year'] >= 2020) & (dois['year'] <= 2023)]

    output_path = DATA_DIR / "dois_for_shi_evans.csv"
    dois.to_csv(output_path, index=False)

    print(f"  Exported {len(dois):,} unique DOIs to {output_path}")

    return dois


def main():
    """Run DOI integration pipeline."""
    print("=" * 60)
    print("DOI INTEGRATION PIPELINE")
    print("=" * 60)

    # Ensure output directory exists
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Load DOI crosswalk
    df = load_doi_crosswalk()

    # Export DOIs for Shi & Evans
    dois_export = export_dois_for_shi_evans(df)

    # Aggregate by gene-semester
    gene_sem_dois = aggregate_dois_by_gene_semester(df)

    # Save
    output_path = DATA_DIR / "gene_semester_dois.parquet"
    gene_sem_dois.to_parquet(output_path, index=False)
    print(f"\nSaved: {output_path}")

    # Summary
    print("\n" + "-" * 40)
    print("SUMMARY")
    print("-" * 40)
    print(f"Total DOIs (2020-2023): {len(dois_export):,}")
    print(f"Gene-semester records: {len(gene_sem_dois):,}")
    print(f"Unique genes with DOIs: {gene_sem_dois['gene_id'].nunique():,}")

    print("\n" + "=" * 60)
    print("DOI INTEGRATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
