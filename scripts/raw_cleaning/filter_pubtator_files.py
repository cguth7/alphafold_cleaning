#!/usr/bin/env python3
"""
Filter PubTator files to keep only essential columns.
- gene2pubtator3.txt: keep pmid (col 1) and gene_id (col 3)  
- disease2pubtator3.txt: keep pmid (col 1) and mesh_id (col 3)
"""

import os
import sys

def filter_gene_file(input_path, output_path):
    """Filter gene PubTator file to keep pmid and gene_id only."""
    print(f"Filtering gene file: {input_path}")
    print(f"Output: {output_path}")
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        # Write header
        outfile.write("pmid\tgene_id\n")
        
        line_count = 0
        for line in infile:
            line_count += 1
            if line_count % 10000000 == 0:
                print(f"Processed {line_count:,} lines...")
                
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pmid = parts[0]
                gene_id = parts[2]
                outfile.write(f"{pmid}\t{gene_id}\n")
    
    print(f"Gene file processing complete. Total lines: {line_count:,}")

def filter_disease_file(input_path, output_path):
    """Filter disease PubTator file to keep pmid and mesh_id only."""
    print(f"Filtering disease file: {input_path}")
    print(f"Output: {output_path}")
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        # Write header
        outfile.write("pmid\tmesh_id\n")
        
        line_count = 0
        for line in infile:
            line_count += 1
            if line_count % 10000000 == 0:
                print(f"Processed {line_count:,} lines...")
                
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pmid = parts[0]
                mesh_id = parts[2]
                outfile.write(f"{pmid}\t{mesh_id}\n")
    
    print(f"Disease file processing complete. Total lines: {line_count:,}")

def main():
    # Define paths
    raw_data_dir = "/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/raw_data"
    cleaned_dir = "/Users/maxguthmann/Downloads/Development/Work/Alphafold_2/cleaned"
    
    gene_input = os.path.join(raw_data_dir, "gene2pubtator3.txt")
    gene_output = os.path.join(cleaned_dir, "gene_pmid_id.txt")
    
    disease_input = os.path.join(raw_data_dir, "disease2pubtator3.txt") 
    disease_output = os.path.join(cleaned_dir, "disease_pmid_mesh.txt")
    
    # Check if input files exist
    if not os.path.exists(gene_input):
        print(f"Error: Gene input file not found: {gene_input}")
        sys.exit(1)
        
    if not os.path.exists(disease_input):
        print(f"Error: Disease input file not found: {disease_input}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(cleaned_dir, exist_ok=True)
    
    # Process files
    print("=== Starting PubTator file filtering ===")
    
    print("\n1. Processing gene file...")
    filter_gene_file(gene_input, gene_output)
    
    print("\n2. Processing disease file...")
    filter_disease_file(disease_input, disease_output)
    
    print("\n=== File filtering complete! ===")
    
    # Report file sizes
    print("\nFile size comparison:")
    for filename, path in [
        ("Original gene file", gene_input),
        ("Filtered gene file", gene_output), 
        ("Original disease file", disease_input),
        ("Filtered disease file", disease_output)
    ]:
        if os.path.exists(path):
            size_mb = os.path.getsize(path) / (1024 * 1024)
            print(f"{filename}: {size_mb:.1f} MB")

if __name__ == "__main__":
    main()