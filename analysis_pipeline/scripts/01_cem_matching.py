"""
01_cem_matching.py - Coarsened Exact Matching for AlphaFold2 Impact Analysis

This script replicates Danilo's CEM matching approach in Python.

Treatment definition:
    treated = 1 if num_deposits == 0 (genes with NO prior PDB structures)

CEM matching variables (pre-treatment, strictly before July 2021):
    - pre_mean_papers: mean publications in pre-period
    - pre_mean_newcom: mean newcomer papers in pre-period
    - pre_mean_veteran: mean veteran papers in pre-period
    - pre_mean_top10: mean top-10% papers in pre-period
    - pre_sd_papers: std dev of papers in pre-period

NOTE: pLDDT is NOT used for matching (per Matteo's instruction - it's an AF2-generated variable)

Output:
    - data/matched_panel.parquet: Full panel with matched genes and CEM weights
    - data/gene_level_features.parquet: Gene-level pre-treatment features
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Configuration
TREATMENT_YM = 738  # July 2021 in Stata ym format
N_QUANTILE_BINS = 10  # Number of bins for CEM coarsening

# Paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "analysis_pipeline" / "data"
FINAL_PANEL = BASE_DIR / "final" / "final_panel_CLEAN.dta"


def load_panel():
    """Load the final panel data."""
    print(f"Loading data from {FINAL_PANEL}")
    df = pd.read_stata(FINAL_PANEL)
    print(f"  Loaded {len(df):,} observations, {df['gene_id'].nunique():,} genes")
    return df


def create_time_variables(df):
    """Create time-related variables matching Danilo's approach."""
    df = df.copy()

    # Sequential month index
    ym_min = df['ym'].min()
    df['ym_seq'] = df['ym'] - ym_min + 1

    # Treatment sequence
    treatment_seq = TREATMENT_YM - ym_min + 1

    # Treatment group: genes with NO prior deposits
    df['treated'] = (df['num_deposits'] == 0).astype(int)

    # Post-treatment indicator
    df['post'] = (df['ym_seq'] >= treatment_seq).astype(int)

    # Quarter and semester indices
    df['quarter'] = ((df['ym_seq'] - 1) // 3).astype(int)
    df['semester'] = ((df['ym_seq'] - 1) // 6).astype(int)

    # Relative time (event time)
    treatment_quarter = (treatment_seq - 1) // 3
    treatment_semester = (treatment_seq - 1) // 6

    df['rel_quarter'] = df['quarter'] - treatment_quarter
    df['rel_semester'] = df['semester'] - treatment_semester

    print(f"  Treatment month (ym_seq): {treatment_seq}")
    print(f"  Treatment quarter: {treatment_quarter}, semester: {treatment_semester}")
    print(f"  Treated genes: {df.groupby('gene_id')['treated'].first().sum():,}")
    print(f"  Control genes: {(~df.groupby('gene_id')['treated'].first().astype(bool)).sum():,}")

    return df, treatment_seq


def compute_pre_treatment_features(df, treatment_seq):
    """Compute gene-level pre-treatment summary statistics."""
    print("Computing pre-treatment features...")

    # Filter to strictly pre-treatment period
    pre_df = df[df['ym_seq'] < treatment_seq].copy()

    # Compute gene-level summaries
    gene_features = pre_df.groupby('gene_id').agg({
        'n_papers': ['mean', 'std'],
        'n_newcomer_papers': 'mean',
        'n_veteran_papers': 'mean',
        'n_top10_y': 'mean',
    }).reset_index()

    # Flatten column names
    gene_features.columns = [
        'gene_id',
        'pre_mean_papers', 'pre_sd_papers',
        'pre_mean_newcom', 'pre_mean_veteran', 'pre_mean_top10'
    ]

    # Fill NaN standard deviations with 0
    gene_features['pre_sd_papers'] = gene_features['pre_sd_papers'].fillna(0)

    # Get gene-level treatment status and pLDDT (for reference, not matching)
    gene_meta = df.groupby('gene_id').agg({
        'treated': 'first',
        'average_plddt': 'first',
        'gene_name': 'first',
        'protein_id': 'first',
    }).reset_index()

    # Merge
    gene_features = gene_features.merge(gene_meta, on='gene_id', how='left')

    # Fill any missing features with 0
    for col in ['pre_mean_papers', 'pre_sd_papers', 'pre_mean_newcom',
                'pre_mean_veteran', 'pre_mean_top10']:
        gene_features[col] = gene_features[col].fillna(0)

    gene_features['average_plddt'] = gene_features['average_plddt'].fillna(0)

    print(f"  Computed features for {len(gene_features):,} genes")

    return gene_features


def coarsened_exact_matching(gene_features, n_bins=N_QUANTILE_BINS):
    """
    Implement Coarsened Exact Matching (CEM).

    CEM works by:
    1. Coarsening continuous covariates into bins
    2. Creating strata from the Cartesian product of bins
    3. Keeping only strata with both treated and control units
    4. Computing weights to balance treated/control within strata
    """
    print(f"Running CEM with {n_bins} quantile bins...")

    df = gene_features.copy()

    # Variables to match on (NOT including pLDDT per Matteo's instruction)
    match_vars = ['pre_mean_papers', 'pre_mean_newcom', 'pre_mean_veteran',
                  'pre_mean_top10', 'pre_sd_papers']

    # Coarsen each variable into quantile bins
    for var in match_vars:
        # Use qcut with duplicates='drop' to handle tied values
        try:
            df[f'cem_{var}'] = pd.qcut(df[var], q=n_bins, labels=False, duplicates='drop')
        except ValueError:
            # If too few unique values, use all unique values as bins
            df[f'cem_{var}'] = pd.cut(df[var], bins=min(n_bins, df[var].nunique()),
                                       labels=False, duplicates='drop')
        df[f'cem_{var}'] = df[f'cem_{var}'].fillna(0).astype(int)

    # Create strata as combination of all bins
    cem_cols = [f'cem_{v}' for v in match_vars]
    df['stratum'] = df[cem_cols].astype(str).agg('_'.join, axis=1)

    # Count treated and control in each stratum
    stratum_counts = df.groupby('stratum').agg({
        'treated': ['sum', 'count']
    }).reset_index()
    stratum_counts.columns = ['stratum', 'n_treated', 'n_total']
    stratum_counts['n_control'] = stratum_counts['n_total'] - stratum_counts['n_treated']

    # Keep only strata with both treated AND control (common support)
    common_support = stratum_counts[
        (stratum_counts['n_treated'] > 0) & (stratum_counts['n_control'] > 0)
    ]['stratum'].tolist()

    print(f"  Total strata: {len(stratum_counts):,}")
    print(f"  Strata with common support: {len(common_support):,}")

    # Filter to common support
    df_matched = df[df['stratum'].isin(common_support)].copy()

    # Merge stratum counts
    df_matched = df_matched.merge(
        stratum_counts[['stratum', 'n_treated', 'n_control']],
        on='stratum', how='left'
    )

    # Compute CEM weights (ATT-style)
    # Treated units get weight 1
    # Control units get weight (n_treated / n_control) within stratum
    df_matched['cem_weight'] = np.where(
        df_matched['treated'] == 1,
        1.0,
        df_matched['n_treated'] / df_matched['n_control']
    )

    # Normalize so total treated weight = total control weight
    sum_treated = df_matched.loc[df_matched['treated'] == 1, 'cem_weight'].sum()
    sum_control = df_matched.loc[df_matched['treated'] == 0, 'cem_weight'].sum()

    df_matched.loc[df_matched['treated'] == 0, 'cem_weight'] *= (sum_treated / sum_control)

    n_treated = (df_matched['treated'] == 1).sum()
    n_control = (df_matched['treated'] == 0).sum()

    print(f"  Matched genes: {len(df_matched):,}")
    print(f"    Treated: {n_treated:,}")
    print(f"    Control: {n_control:,}")

    return df_matched


def create_matched_panel(panel_df, matched_genes):
    """Merge matched genes back to full panel."""
    print("Creating matched panel...")

    # Keep only matched genes
    matched_gene_ids = matched_genes['gene_id'].unique()

    panel_matched = panel_df[panel_df['gene_id'].isin(matched_gene_ids)].copy()

    # Merge CEM weights
    weight_cols = ['gene_id', 'cem_weight', 'stratum', 'treated']
    panel_matched = panel_matched.drop(columns=['treated'], errors='ignore')
    panel_matched = panel_matched.merge(
        matched_genes[weight_cols],
        on='gene_id',
        how='left'
    )

    print(f"  Matched panel: {len(panel_matched):,} observations")
    print(f"  Genes: {panel_matched['gene_id'].nunique():,}")

    return panel_matched


def main():
    """Run the CEM matching pipeline."""
    print("=" * 60)
    print("CEM MATCHING PIPELINE")
    print("=" * 60)

    # Ensure output directory exists
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    df = load_panel()

    # Create time variables
    df, treatment_seq = create_time_variables(df)

    # Compute pre-treatment features
    gene_features = compute_pre_treatment_features(df, treatment_seq)

    # Run CEM
    matched_genes = coarsened_exact_matching(gene_features)

    # Create matched panel
    panel_matched = create_matched_panel(df, matched_genes)

    # Save outputs
    print("\nSaving outputs...")

    gene_features.to_parquet(DATA_DIR / "gene_level_features.parquet", index=False)
    print(f"  Saved: {DATA_DIR / 'gene_level_features.parquet'}")

    matched_genes.to_parquet(DATA_DIR / "matched_genes.parquet", index=False)
    print(f"  Saved: {DATA_DIR / 'matched_genes.parquet'}")

    panel_matched.to_parquet(DATA_DIR / "matched_panel_monthly.parquet", index=False)
    print(f"  Saved: {DATA_DIR / 'matched_panel_monthly.parquet'}")

    print("\n" + "=" * 60)
    print("CEM MATCHING COMPLETE")
    print("=" * 60)

    return panel_matched, matched_genes


if __name__ == "__main__":
    main()
