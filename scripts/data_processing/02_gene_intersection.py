#!/usr/bin/env python3
"""
Phase 2: Gene-Disease Intersection

Filter gene PubTator data to keep only:
1. Gene mentions in disease-relevant publications (from Phase 1)  
2. Genes that exist in the master protein dataset

Input:  
- cleaned/gene_pmid_id.txt (72.9M rows, 1.0GB)
- processed/disease_relevant_pmids.parquet (21.8M PMIDs)
- intermediate_data/protein_data_master_dta.dta (1.2M genes)

Output: processed/gene_disease_intersect.parquet

Processing steps:
1. Load cleaned gene data (pmid, gene_id)
2. Inner merge with disease-relevant PMIDs 
3. Filter to master gene list
4. Document merge success rates and excluded entries
5. Save with full statistics
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def load_gene_data(file_path):
    """Load cleaned gene PubTator data."""
    print(f"Loading gene data from: {file_path}")
    
    # Read tab-separated file as strings first to handle semicolon-separated IDs
    df = pd.read_csv(file_path, sep='\t', dtype={'pmid': 'int32', 'gene_id': 'str'})
    
    initial_count = len(df)
    print(f"  Loaded {initial_count:,} gene mentions")
    
    # Check for semicolon-separated gene IDs
    semicolon_mask = df['gene_id'].str.contains(';', na=False)
    semicolon_count = semicolon_mask.sum()
    
    if semicolon_count > 0:
        print(f"  Found {semicolon_count:,} entries with semicolon-separated gene IDs ({semicolon_count/initial_count*100:.2f}%)")
        print(f"  Sample problematic entries: {df[semicolon_mask]['gene_id'].head(3).tolist()}")
        
        # Drop semicolon-separated entries (as per your Stata comment)
        df = df[~semicolon_mask].copy()
        print(f"  Dropped semicolon entries, remaining: {len(df):,} ({len(df)/initial_count*100:.2f}%)")
    
    # Convert gene_id to int64 to handle large gene IDs
    df['gene_id'] = df['gene_id'].astype('int64')
    
    print(f"  Unique PMIDs: {df['pmid'].nunique():,}")
    print(f"  Unique genes: {df['gene_id'].nunique():,}")
    print(f"  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    return df

def load_disease_pmids(file_path):
    """Load disease-relevant PMIDs from Phase 1."""
    print(f"Loading disease-relevant PMIDs from: {file_path}")
    
    df = pd.read_parquet(file_path)
    pmid_set = set(df['pmid'].values)
    
    print(f"  Loaded {len(pmid_set):,} disease-relevant PMIDs")
    
    return pmid_set

def intersect_gene_disease(df_gene, disease_pmid_set):
    """Inner merge gene data with disease-relevant PMIDs."""
    print("Intersecting gene mentions with disease-relevant publications...")
    
    initial_count = len(df_gene)
    initial_pmids = df_gene['pmid'].nunique()
    initial_genes = df_gene['gene_id'].nunique()
    
    # Filter to disease-relevant PMIDs
    df_filtered = df_gene[df_gene['pmid'].isin(disease_pmid_set)].copy()
    
    final_count = len(df_filtered)
    final_pmids = df_filtered['pmid'].nunique()
    final_genes = df_filtered['gene_id'].nunique()
    
    print(f"  Initial gene mentions: {initial_count:,}")
    print(f"  Initial PMIDs: {initial_pmids:,}")
    print(f"  Initial genes: {initial_genes:,}")
    print(f"  Final gene mentions: {final_count:,} ({final_count/initial_count*100:.2f}% retained)")
    print(f"  Final PMIDs: {final_pmids:,} ({final_pmids/initial_pmids*100:.2f}% retained)")
    print(f"  Final genes: {final_genes:,} ({final_genes/initial_genes*100:.2f}% retained)")
    print(f"  Avg gene mentions per paper: {final_count/final_pmids:.2f}")
    
    return df_filtered, {
        'initial_mentions': initial_count,
        'final_mentions': final_count,
        'initial_pmids': initial_pmids,
        'final_pmids': final_pmids,
        'initial_genes': initial_genes,
        'final_genes': final_genes,
        'mention_retention': final_count/initial_count,
        'pmid_retention': final_pmids/initial_pmids,
        'gene_retention': final_genes/initial_genes
    }

def load_master_genes(file_path):
    """Load master gene list from protein data."""
    print(f"Loading master gene list from: {file_path}")
    
    df = pd.read_stata(file_path)
    
    # Extract unique gene IDs  
    master_genes = df['geneid'].dropna().astype('int64').unique()
    master_gene_set = set(master_genes)
    
    print(f"  Loaded {len(master_gene_set):,} master genes")
    
    # Sample of master genes for documentation
    sample_genes = sorted(master_genes)[:10]
    print(f"  Sample master genes: {sample_genes}")
    
    return master_gene_set, {
        'total_master_genes': len(master_gene_set),
        'sample_genes': sample_genes
    }

def filter_to_master_genes(df, master_gene_set):
    """Filter gene-disease intersected data to master gene list."""
    print("Filtering to master gene list...")
    
    initial_count = len(df)
    initial_genes = df['gene_id'].nunique()
    
    # Check which genes are in master list
    df['in_master'] = df['gene_id'].isin(master_gene_set)
    
    # Get statistics before filtering
    genes_in_master = df['gene_id'][df['in_master']].nunique()
    genes_not_in_master = df['gene_id'][~df['in_master']].nunique()
    
    # Filter to master genes only
    df_filtered = df[df['in_master']].copy().drop('in_master', axis=1)
    
    final_count = len(df_filtered)
    final_genes = df_filtered['gene_id'].nunique()
    
    print(f"  Initial mentions: {initial_count:,}")
    print(f"  Initial genes: {initial_genes:,}")
    print(f"  Genes in master list: {genes_in_master:,}")
    print(f"  Genes NOT in master list: {genes_not_in_master:,}")
    print(f"  Final mentions: {final_count:,} ({final_count/initial_count*100:.2f}% retained)")
    print(f"  Final genes: {final_genes:,} ({final_genes/initial_genes*100:.2f}% retained)")
    
    # Sample of excluded genes for documentation
    excluded_genes = df['gene_id'][~df['in_master']].unique()[:10]
    print(f"  Sample excluded genes: {excluded_genes.tolist()}")
    
    # Critical check: did we lose a lot of data?
    retention_rate = final_count / initial_count
    if retention_rate < 0.8:
        print(f"  ⚠️  WARNING: Low retention rate ({retention_rate*100:.1f}%)")
        print(f"     This suggests many genes in PubTator are not in master list")
    else:
        print(f"  ✅ Good retention rate ({retention_rate*100:.1f}%)")
    
    return df_filtered, {
        'initial_mentions': initial_count,
        'final_mentions': final_count,
        'initial_genes': initial_genes,
        'final_genes': final_genes,
        'genes_in_master': genes_in_master,
        'genes_not_in_master': genes_not_in_master,
        'mention_retention': retention_rate,
        'gene_retention': final_genes/initial_genes,
        'excluded_sample': excluded_genes[:5].tolist()
    }

def save_results(df, output_path, intersection_stats, master_stats):
    """Save results with comprehensive metadata."""
    print(f"Saving results to: {output_path}")
    
    # Save main output
    df.to_parquet(output_path, index=False)
    
    file_size = os.path.getsize(output_path) / 1024**2
    print(f"  Saved {len(df):,} gene-disease intersection records")
    print(f"  File size: {file_size:.1f} MB")
    
    # Save processing statistics
    stats_path = output_path.parent / 'phase2_processing_stats.txt'
    with open(stats_path, 'w') as f:
        f.write("Phase 2: Gene-Disease Intersection Processing Statistics\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Input files:\n")
        f.write("  - cleaned/gene_pmid_id.txt (Gene PubTator data)\n")
        f.write("  - processed/disease_relevant_pmids.parquet (Phase 1 output)\n")
        f.write("  - intermediate_data/protein_data_master_dta.dta (Master genes)\n\n")
        
        f.write(f"Output file: {output_path.name}\n\n")
        
        f.write("STEP 1: Gene-Disease Publication Intersection\n")
        f.write("-" * 45 + "\n")
        f.write(f"Initial gene mentions: {intersection_stats['initial_mentions']:,}\n")
        f.write(f"Disease-relevant mentions: {intersection_stats['final_mentions']:,}\n")
        f.write(f"Mention retention: {intersection_stats['mention_retention']*100:.2f}%\n\n")
        
        f.write(f"Initial unique PMIDs: {intersection_stats['initial_pmids']:,}\n")
        f.write(f"Disease-relevant PMIDs: {intersection_stats['final_pmids']:,}\n")
        f.write(f"PMID retention: {intersection_stats['pmid_retention']*100:.2f}%\n\n")
        
        f.write(f"Initial unique genes: {intersection_stats['initial_genes']:,}\n")
        f.write(f"Disease-relevant genes: {intersection_stats['final_genes']:,}\n")
        f.write(f"Gene retention: {intersection_stats['gene_retention']*100:.2f}%\n\n")
        
        f.write("STEP 2: Master Gene List Filtering\n")
        f.write("-" * 35 + "\n")
        f.write(f"Master genes available: {master_stats['total_master_genes']:,}\n")
        f.write(f"Genes in master list: {master_stats['genes_in_master']:,}\n")
        f.write(f"Genes NOT in master: {master_stats['genes_not_in_master']:,}\n")
        f.write(f"Final gene retention: {master_stats['gene_retention']*100:.2f}%\n")
        f.write(f"Final mention retention: {master_stats['mention_retention']*100:.2f}%\n\n")
        
        f.write("Final Output Summary\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Unique PMIDs: {df['pmid'].nunique():,}\n")
        f.write(f"Unique genes: {df['gene_id'].nunique():,}\n")
        f.write(f"Avg mentions per paper: {len(df)/df['pmid'].nunique():.2f}\n")
        f.write(f"Output file size: {file_size:.1f} MB\n\n")
        
        f.write("Sample excluded genes (not in master):\n")
        for i, gene in enumerate(master_stats['excluded_sample'], 1):
            f.write(f"  {i}. Gene ID: {gene}\n")
        
        if master_stats['mention_retention'] < 0.8:
            f.write(f"\n⚠️  WARNING: Low retention rate ({master_stats['mention_retention']*100:.1f}%)\n")
            f.write("   Many genes in PubTator are missing from master protein dataset.\n")
            f.write("   Consider investigating master gene list completeness.\n")
    
    print(f"  Processing stats saved to: {stats_path}")

def main():
    # Define paths
    base_dir = Path("/Users/maxguthmann/Downloads/Development/Work/Alphafold_2")
    
    gene_input = base_dir / "cleaned" / "gene_pmid_id.txt"
    disease_pmids = base_dir / "processed" / "disease_relevant_pmids.parquet"
    master_genes = base_dir / "intermediate_data" / "protein_data_master_dta.dta"
    output_file = base_dir / "processed" / "gene_disease_intersect.parquet"
    
    print("=== Phase 2: Gene-Disease Intersection ===\n")
    
    # Step 1: Load gene data
    df_gene = load_gene_data(gene_input)
    
    # Step 2: Load disease-relevant PMIDs
    disease_pmid_set = load_disease_pmids(disease_pmids)
    
    # Step 3: Intersect gene mentions with disease publications
    df_intersect, intersection_stats = intersect_gene_disease(df_gene, disease_pmid_set)
    
    # Step 4: Load master gene list
    master_gene_set, master_info = load_master_genes(master_genes)
    
    # Step 5: Filter to master genes
    df_filtered, master_filter_stats = filter_to_master_genes(df_intersect, master_gene_set)
    
    # Combine master stats
    master_stats = {**master_info, **master_filter_stats}
    
    # Step 6: Save results
    save_results(df_filtered, output_file, intersection_stats, master_stats)
    
    print(f"\n=== Phase 2 Complete ===")
    print(f"Output: {output_file}")
    print(f"Gene-disease intersection records: {len(df_filtered):,}")
    print(f"Unique genes: {df_filtered['gene_id'].nunique():,}")
    print(f"Unique PMIDs: {df_filtered['pmid'].nunique():,}")

if __name__ == "__main__":
    main()