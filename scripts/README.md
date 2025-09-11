# Scripts Directory

This directory contains all scripts for data processing and analysis, organized by processing stage.

## Directory Structure

```
scripts/
├── raw_cleaning/          # Scripts for initial data cleaning from raw_data/
├── data_processing/       # Scripts for multi-phase data processing pipeline
└── README.md             # This file
```

## Usage

All Python scripts should be run using the system Python interpreter:
```bash
/usr/bin/python3 script_name.py
```

## Scripts

### raw_cleaning/

#### filter_pubtator_files.py
**Purpose:** Filter PubTator files to keep only essential columns
- Extracts pmid + gene_id from gene2pubtator3.txt  
- Extracts pmid + mesh_id from disease2pubtator3.txt
- Outputs cleaned files to cleaned/ directory

**Usage:**
```bash
cd /Users/maxguthmann/Downloads/Development/Work/Alphafold_2
/usr/bin/python3 scripts/raw_cleaning/filter_pubtator_files.py
```

**Input:** raw_data/gene2pubtator3.txt, raw_data/disease2pubtator3.txt  
**Output:** cleaned/gene_pmid_id.txt, cleaned/disease_pmid_mesh.txt

### data_processing/

#### 01_disease_filter.py
**Purpose:** Phase 1 - Filter disease PubTator data for relevant publications  
- Keeps only MESH depth 3/4 diseases OR OMIM diseases
- Creates unique PMID list of disease-relevant publications
- Comprehensive processing statistics and validation

**Usage:**
```bash
cd /Users/maxguthmann/Downloads/Development/Work/Alphafold_2
/usr/bin/python3 scripts/data_processing/01_disease_filter.py
```

**Input:** cleaned/disease_pmid_mesh.txt  
**Output:** processed/disease_relevant_pmids.parquet, phase1_processing_stats.txt

#### 02_gene_intersection.py
**Purpose:** Phase 2 - Create gene-disease intersection dataset  
- Inner merge gene PubTator data with disease-relevant PMIDs from Phase 1
- Filter to master protein gene list 
- Handle semicolon-separated gene IDs and large gene ID values
- Comprehensive merge statistics and data quality checks

**Usage:**
```bash
cd /Users/maxguthmann/Downloads/Development/Work/Alphafold_2  
/usr/bin/python3 scripts/data_processing/02_gene_intersection.py
```

**Input:** cleaned/gene_pmid_id.txt, processed/disease_relevant_pmids.parquet, intermediate_data/protein_data_master_dta.dta  
**Output:** processed/gene_disease_intersect.parquet, phase2_processing_stats.txt