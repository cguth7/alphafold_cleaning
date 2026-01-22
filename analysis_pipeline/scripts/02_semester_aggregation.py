"""
02_semester_aggregation.py - Create Semester-Level Panel with Delta Variables

This script:
1. Aggregates monthly data to semester level
2. Creates ΔY (first difference) variables
3. Creates asinh transformations and their differences

Output:
    - data/matched_panel_semester.parquet: Semester-level panel ready for event studies
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "analysis_pipeline" / "data"

# Outcomes to process
OUTCOMES = ['n_papers', 'n_newcomer_papers', 'n_veteran_papers', 'n_top10_y']


def load_matched_panel():
    """Load the matched monthly panel."""
    path = DATA_DIR / "matched_panel_monthly.parquet"
    print(f"Loading matched panel from {path}")
    df = pd.read_parquet(path)
    print(f"  Loaded {len(df):,} observations")
    return df


def aggregate_to_semester(df):
    """Aggregate monthly data to semester level."""
    print("Aggregating to semester level...")

    # Group by gene and semester
    agg_dict = {
        # Sum outcome variables
        'n_papers': 'sum',
        'n_newcomer_papers': 'sum',
        'n_veteran_papers': 'sum',
        'n_top10_y': 'sum',
        'n_top25_y': 'sum',
        'n_top05_y': 'sum',
        # Keep first for gene-level variables
        'treated': 'first',
        'cem_weight': 'first',
        'stratum': 'first',
        'gene_name': 'first',
        'protein_id': 'first',
        'average_plddt': 'first',
        'num_deposits': 'first',
        # Keep relative time
        'rel_semester': 'first',
        'post': 'max',  # 1 if any month is post
        # Disease metrics
        'unique_mesh_count': 'sum',
        'new_mesh_count': 'sum',
    }

    # Only aggregate columns that exist
    agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}

    df_sem = df.groupby(['gene_id', 'semester']).agg(agg_dict).reset_index()

    print(f"  Created {len(df_sem):,} gene-semester observations")
    print(f"  Semesters: {df_sem['semester'].min()} to {df_sem['semester'].max()}")
    print(f"  Relative semesters: {df_sem['rel_semester'].min()} to {df_sem['rel_semester'].max()}")

    return df_sem


def create_delta_variables(df):
    """Create first-difference (ΔY) and asinh variables."""
    print("Creating delta and asinh variables...")

    df = df.sort_values(['gene_id', 'semester']).copy()

    for y in OUTCOMES:
        if y not in df.columns:
            print(f"  Skipping {y} (not in data)")
            continue

        # Asinh transformation (handles zeros better than log)
        df[f'asinh_{y}'] = np.arcsinh(df[y])

        # First differences within gene
        df[f'D_{y}'] = df.groupby('gene_id')[y].diff()
        df[f'D_asinh_{y}'] = df.groupby('gene_id')[f'asinh_{y}'].diff()

    # Count non-null deltas
    n_valid = df['D_n_papers'].notna().sum()
    print(f"  Valid delta observations: {n_valid:,}")

    return df


def create_event_study_variables(df):
    """Create variables needed for event study regressions."""
    print("Creating event study variables...")

    # Factor-variable safe relative time (non-negative)
    min_sem = df['rel_semester'].min()
    df['rel_sem_fv'] = df['rel_semester'] - min_sem

    # Linear time trend
    df['t_sem'] = df['rel_semester'].astype(float)

    # Pre-period indicator
    df['pre'] = (df['rel_semester'] < 0).astype(int)

    print(f"  rel_sem_fv range: {df['rel_sem_fv'].min()} to {df['rel_sem_fv'].max()}")

    return df


def main():
    """Run semester aggregation pipeline."""
    print("=" * 60)
    print("SEMESTER AGGREGATION PIPELINE")
    print("=" * 60)

    # Load matched panel
    df = load_matched_panel()

    # Aggregate to semester
    df_sem = aggregate_to_semester(df)

    # Create delta variables
    df_sem = create_delta_variables(df_sem)

    # Create event study variables
    df_sem = create_event_study_variables(df_sem)

    # Save
    output_path = DATA_DIR / "matched_panel_semester.parquet"
    df_sem.to_parquet(output_path, index=False)
    print(f"\nSaved: {output_path}")

    # Summary statistics
    print("\n" + "-" * 40)
    print("SUMMARY STATISTICS")
    print("-" * 40)

    for y in OUTCOMES:
        if y in df_sem.columns:
            print(f"\n{y}:")
            print(f"  Mean: {df_sem[y].mean():.3f}")
            print(f"  Std:  {df_sem[y].std():.3f}")
            if f'D_{y}' in df_sem.columns:
                print(f"  Δ Mean: {df_sem[f'D_{y}'].mean():.4f}")
                print(f"  Δ Std:  {df_sem[f'D_{y}'].std():.4f}")

    print("\n" + "=" * 60)
    print("SEMESTER AGGREGATION COMPLETE")
    print("=" * 60)

    return df_sem


if __name__ == "__main__":
    main()
