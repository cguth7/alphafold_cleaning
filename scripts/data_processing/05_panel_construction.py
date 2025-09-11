#!/usr/bin/env python3
"""
Phase 5: Panel Construction

Create the final balanced Gene-Month panel for AlphaFold impact analysis.

Following the Stata approach:
1. Aggregate citation-enriched data to Gene-Month level with all metrics
2. Balance the panel across all gene-month combinations (2020-2023)
3. Merge in unique disease metrics from intermediate data
4. Add master gene metadata (protein info, gene families, etc.)
5. Fill missing values with appropriate defaults (0 for counts)

Input:  
- processed/citation_enriched.parquet (13.7M records)
- intermediate_data/gene_month_unique_mesh_2015.dta
- intermediate_data/gene_month_new_unique_diseases_2015_no_level_5.dta  
- intermediate_data/protein_data_master_dta.dta

Output: final/final_balanced_panel.parquet

Processing steps match your Stata workflow for perfect replication.
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from itertools import product

def load_citation_enriched_data(file_path):
    """Load citation-enriched data from Phase 4."""
    print(f"Loading citation-enriched data from: {file_path}")
    
    df = pd.read_parquet(file_path)
    
    print(f"  Loaded {len(df):,} citation-enriched records")
    print(f"  Time span: {df['year'].min()}-{df['year'].max()}")
    print(f"  Unique PMIDs: {df['pmid'].nunique():,}")
    print(f"  Unique genes: {df['gene_id'].nunique():,}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def aggregate_to_gene_month(df):
    """Aggregate to Gene-Month level following Stata approach."""
    print("Aggregating to Gene-Month level...")
    
    # Step 1: Handle multiple gene mentions per paper (like Stata bysort Gene_ID ym pmid: keep if _n == 1)
    print("  Step 1: One row per Gene-Month-PMID combination...")
    
    # For each gene-month-pmid, take max of any flags (handles multiple author records per paper)
    paper_gene_level = (df.groupby(['gene_id', 'ym', 'pmid'])
                       .agg({
                           'year': 'first',
                           'month': 'first',
                           'quarter': 'first', 
                           'year_quarter': 'first',
                           'bimonth': 'first',
                           'year_bimonth': 'first',
                           'icite_cites': 'first',
                           # High-citation flags (max ensures if any author hit threshold, paper counts)
                           'hi25_year': 'max', 'hi10_year': 'max', 'hi05_year': 'max',
                           'hi25_quarter': 'max', 'hi10_quarter': 'max', 'hi05_quarter': 'max', 
                           'hi25_bimonth': 'max', 'hi10_bimonth': 'max', 'hi05_bimonth': 'max',
                           # Author novelty flags (max means paper has at least one newcomer/veteran)
                           'newcomer_author': 'max',
                           'not_newcomer_author': 'max'
                       })
                       .reset_index())
    
    print(f"    Unique gene-month-paper combinations: {len(paper_gene_level):,}")
    
    # Step 2: Collapse to Gene-Month level (like Stata collapse command)
    print("  Step 2: Collapsing to Gene-Month with counts and sums...")
    
    gene_month_panel = (paper_gene_level.groupby(['gene_id', 'ym'])
                       .agg({
                           'year': 'first',
                           'month': 'first', 
                           'quarter': 'first',
                           'year_quarter': 'first',
                           'bimonth': 'first',
                           'year_bimonth': 'first',
                           # Paper counts
                           'pmid': 'count',  # n_papers
                           # High-citation paper counts (annual)
                           'hi25_year': 'sum',   # n_top25_y  
                           'hi10_year': 'sum',   # n_top10_y
                           'hi05_year': 'sum',   # n_top05_y
                           # High-citation paper counts (quarterly)
                           'hi25_quarter': 'sum',   # n_top25_q
                           'hi10_quarter': 'sum',   # n_top10_q  
                           'hi05_quarter': 'sum',   # n_top05_q
                           # High-citation paper counts (bi-monthly)
                           'hi25_bimonth': 'sum',   # n_top25_b2
                           'hi10_bimonth': 'sum',   # n_top10_b2
                           'hi05_bimonth': 'sum',   # n_top05_b2
                           # Author novelty counts
                           'newcomer_author': 'sum',       # papers with newcomer authors
                           'not_newcomer_author': 'sum'    # papers with veteran authors
                       })
                       .reset_index())
    
    # Rename columns to match Stata output
    gene_month_panel = gene_month_panel.rename(columns={
        'pmid': 'n_papers',
        'hi25_year': 'n_top25_y', 'hi10_year': 'n_top10_y', 'hi05_year': 'n_top05_y',
        'hi25_quarter': 'n_top25_q', 'hi10_quarter': 'n_top10_q', 'hi05_quarter': 'n_top05_q', 
        'hi25_bimonth': 'n_top25_b2', 'hi10_bimonth': 'n_top10_b2', 'hi05_bimonth': 'n_top05_b2',
        'newcomer_author': 'n_newcomer_papers',
        'not_newcomer_author': 'n_veteran_papers'
    })
    
    print(f"  Final gene-month observations: {len(gene_month_panel):,}")
    print(f"  Unique genes: {gene_month_panel['gene_id'].nunique():,}")
    print(f"  Time span: {gene_month_panel['ym'].min()}-{gene_month_panel['ym'].max()}")
    
    # Summary statistics
    print(f"  Total papers: {gene_month_panel['n_papers'].sum():,}")
    print(f"  Mean papers per gene-month: {gene_month_panel['n_papers'].mean():.2f}")
    print(f"  Gene-months with activity: {(gene_month_panel['n_papers'] > 0).sum():,}")
    
    return gene_month_panel

def create_balanced_panel(df_aggregated):
    """Create balanced panel with all gene-month combinations."""
    print("Creating balanced panel...")
    
    # Get all unique genes and time periods
    all_genes = sorted(df_aggregated['gene_id'].unique())
    all_ym = sorted(df_aggregated['ym'].unique())
    
    print(f"  Genes: {len(all_genes):,}")
    print(f"  Time periods: {len(all_ym)} months ({min(all_ym)}-{max(all_ym)})")
    print(f"  Expected balanced panel size: {len(all_genes) * len(all_ym):,}")
    
    # Create full cartesian product  
    balanced_index = pd.DataFrame(list(product(all_genes, all_ym)), 
                                columns=['gene_id', 'ym'])
    
    print(f"  Created balanced index: {len(balanced_index):,} rows")
    
    # Merge with aggregated data
    balanced_panel = balanced_index.merge(df_aggregated, on=['gene_id', 'ym'], how='left')
    
    # Fill missing temporal variables
    print("  Filling missing temporal variables...")
    for ym in all_ym:
        mask = balanced_panel['ym'] == ym
        if balanced_panel.loc[mask, 'year'].isna().any():
            # Extract year and month from ym (YYYYMM format)
            year = ym // 100
            month = ym % 100
            
            balanced_panel.loc[mask, 'year'] = year
            balanced_panel.loc[mask, 'month'] = month
            balanced_panel.loc[mask, 'quarter'] = ((month - 1) // 3) + 1
            balanced_panel.loc[mask, 'year_quarter'] = year * 10 + balanced_panel.loc[mask, 'quarter']
            balanced_panel.loc[mask, 'bimonth'] = ((month - 1) // 2) + 1  
            balanced_panel.loc[mask, 'year_bimonth'] = year * 10 + balanced_panel.loc[mask, 'bimonth']
    
    # Fill missing count variables with 0
    count_columns = [
        'n_papers', 'n_top25_y', 'n_top10_y', 'n_top05_y',
        'n_top25_q', 'n_top10_q', 'n_top05_q',
        'n_top25_b2', 'n_top10_b2', 'n_top05_b2',
        'n_newcomer_papers', 'n_veteran_papers'
    ]
    
    for col in count_columns:
        balanced_panel[col] = balanced_panel[col].fillna(0).astype('int32')
    
    print(f"  Balanced panel created: {len(balanced_panel):,} rows")
    print(f"  Non-zero observations: {(balanced_panel['n_papers'] > 0).sum():,} ({(balanced_panel['n_papers'] > 0).mean()*100:.2f}%)")
    
    return balanced_panel

def merge_unique_disease_metrics(balanced_panel, unique_mesh_file, new_unique_file):
    """Merge in unique disease metrics from intermediate data."""
    print("Merging unique disease metrics...")
    
    # Load unique mesh data
    print(f"  Loading unique mesh data: {unique_mesh_file}")
    df_unique = pd.read_stata(unique_mesh_file)
    print(f"    Loaded {len(df_unique):,} unique mesh records")
    print(f"    Columns: {list(df_unique.columns)}")
    
    # Load new unique diseases  
    print(f"  Loading new unique diseases: {new_unique_file}")
    df_new_unique = pd.read_stata(new_unique_file)
    print(f"    Loaded {len(df_new_unique):,} new unique disease records")
    print(f"    Columns: {list(df_new_unique.columns)}")
    
    # Merge unique mesh counts
    print("  Merging unique mesh counts...")
    merged = balanced_panel.merge(df_unique[['Gene_ID', 'ym', 'unique_mesh_count']], 
                                left_on=['gene_id', 'ym'], 
                                right_on=['Gene_ID', 'ym'], 
                                how='left')
    merged = merged.drop('Gene_ID', axis=1)
    merged['unique_mesh_count'] = merged['unique_mesh_count'].fillna(0).astype('int32')
    
    # Merge new unique disease counts  
    print("  Merging new unique disease counts...")
    merged = merged.merge(df_new_unique[['Gene_ID', 'ym', 'new_mesh_count']], 
                        left_on=['gene_id', 'ym'],
                        right_on=['Gene_ID', 'ym'], 
                        how='left')
    merged = merged.drop('Gene_ID', axis=1)
    merged['new_mesh_count'] = merged['new_mesh_count'].fillna(0).astype('int32')
    
    unique_mesh_added = (merged['unique_mesh_count'] > 0).sum()
    new_mesh_added = (merged['new_mesh_count'] > 0).sum()
    
    print(f"  Unique mesh counts added to {unique_mesh_added:,} observations")
    print(f"  New mesh counts added to {new_mesh_added:,} observations")
    
    return merged

def add_master_gene_metadata(balanced_panel, master_file):
    """Add master gene metadata (protein info, etc.)."""
    print("Adding master gene metadata...")
    
    print(f"  Loading master gene data: {master_file}")
    df_master = pd.read_stata(master_file)
    print(f"    Loaded {len(df_master):,} master records")
    
    # Get unique gene metadata
    gene_metadata = (df_master.groupby('geneid')
                    .agg({
                        'protein_id': 'first',
                        'protein_existence': 'first', 
                        'gene_name': 'first',
                        'average_plddt': 'mean'  # Average across time if multiple
                    })
                    .reset_index())
    
    print(f"    Unique genes in master: {len(gene_metadata):,}")
    
    # Merge with balanced panel
    merged = balanced_panel.merge(gene_metadata, 
                                left_on='gene_id', 
                                right_on='geneid', 
                                how='left')
    merged = merged.drop('geneid', axis=1)
    
    genes_with_metadata = merged['protein_id'].notna().sum()
    print(f"  Gene metadata added to {genes_with_metadata:,} observations")
    
    return merged

def save_final_panel(df, output_path):
    """Save final balanced panel with comprehensive metadata."""
    print(f"Saving final balanced panel to: {output_path}")
    
    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save main output
    df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} gene-month observations")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save comprehensive statistics
    stats_path = output_path.parent / 'final_panel_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Final Balanced Gene-Month Panel Statistics\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Panel Structure:\n")
        f.write(f"  Total observations: {len(df):,}\n")
        f.write(f"  Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"  Time periods: {df['ym'].nunique()} months\n")
        f.write(f"  Time span: {df['year'].min()}-{df['year'].max()}\n")
        f.write(f"  Balanced: {'Yes' if len(df) == df['gene_id'].nunique() * df['ym'].nunique() else 'No'}\n\n")
        
        f.write("Activity Summary:\n")
        active_obs = (df['n_papers'] > 0).sum()
        f.write(f"  Observations with publications: {active_obs:,} ({active_obs/len(df)*100:.2f}%)\n")
        f.write(f"  Total publications: {df['n_papers'].sum():,}\n")
        f.write(f"  Mean publications per gene-month: {df['n_papers'].mean():.3f}\n")
        f.write(f"  Max publications in gene-month: {df['n_papers'].max()}\n\n")
        
        f.write("Citation Quality:\n")
        f.write(f"  High-citation papers (top 25% annual): {df['n_top25_y'].sum():,}\n")
        f.write(f"  High-citation papers (top 10% annual): {df['n_top10_y'].sum():,}\n")
        f.write(f"  High-citation papers (top 5% annual): {df['n_top05_y'].sum():,}\n\n")
        
        f.write("Author Dynamics:\n")
        f.write(f"  Papers with newcomer authors: {df['n_newcomer_papers'].sum():,}\n")
        f.write(f"  Papers with veteran authors: {df['n_veteran_papers'].sum():,}\n\n")
        
        f.write("Disease Innovation:\n")
        f.write(f"  Unique mesh diseases: {df['unique_mesh_count'].sum():,}\n")
        f.write(f"  New mesh diseases: {df['new_mesh_count'].sum():,}\n\n")
        
        f.write("Gene Coverage:\n")
        genes_with_metadata = df['protein_id'].notna().sum() // df['ym'].nunique()  
        f.write(f"  Genes with protein metadata: {genes_with_metadata:,}\n")
        f.write(f"  Coverage rate: {genes_with_metadata/df['gene_id'].nunique()*100:.2f}%\n\n")
        
        f.write(f"Output file size: {file_size:.1f} MB\n")
        f.write("Ready for AlphaFold impact analysis!\n")
    
    print(f"  Panel statistics saved to: {stats_path}")
    
    # Quick validation
    expected_size = df['gene_id'].nunique() * df['ym'].nunique()
    if len(df) == expected_size:
        print(f"  ✅ Panel is perfectly balanced ({expected_size:,} observations)")
    else:
        print(f"  ⚠️ Panel size mismatch: {len(df):,} vs expected {expected_size:,}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    citation_input = base_dir / "processed" / "citation_enriched.parquet"
    unique_mesh_file = base_dir / "intermediate_data" / "gene_month_unique_mesh_2015.dta"
    new_unique_file = base_dir / "intermediate_data" / "gene_month_new_unique_diseases_2015_no_level_5.dta"
    master_file = base_dir / "intermediate_data" / "protein_data_master_dta.dta"
    output_file = base_dir / "final" / "final_balanced_panel.parquet"
    
    print("=== Phase 5: Panel Construction ===\n")
    
    # Step 1: Load citation-enriched data
    df_citation = load_citation_enriched_data(citation_input)
    
    # Step 2: Aggregate to gene-month level
    print("\\n" + "="*50)
    df_aggregated = aggregate_to_gene_month(df_citation)
    
    # Step 3: Create balanced panel
    print("\\n" + "="*50)
    balanced_panel = create_balanced_panel(df_aggregated)
    
    # Step 4: Merge unique disease metrics
    print("\\n" + "="*50)
    with_unique_diseases = merge_unique_disease_metrics(balanced_panel, unique_mesh_file, new_unique_file)
    
    # Step 5: Add master gene metadata
    print("\\n" + "="*50)
    final_panel = add_master_gene_metadata(with_unique_diseases, master_file)
    
    # Step 6: Save final panel
    print("\\n" + "="*50)
    save_final_panel(final_panel, output_file)
    
    print(f"\\n=== Phase 5 Complete ===")
    print(f"Final balanced panel: {len(final_panel):,} observations")
    print(f"Ready for AlphaFold impact analysis!")
    print(f"Output: {output_file}")

if __name__ == "__main__":
    main()