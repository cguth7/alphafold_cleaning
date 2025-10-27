# UniProt Citation Data Extraction

**Monthly Publication Analysis for Human Proteins**

This pipeline extracts publication citation data from UniProt for all human proteins (~20,000 entries) and aggregates it by month to analyze research trends over time.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Data Quality](#data-quality)
- [Known Limitations](#known-limitations)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Overview

UniProt provides curated and computationally-mapped publication citations for protein entries. This pipeline:

1. Downloads all human protein entries from UniProt (reviewed + unreviewed)
2. Extracts publication data with dates
3. Aggregates citations by month and year
4. Validates data quality
5. Produces multiple output formats for analysis

**Total expected runtime:** ~1-2 hours for full pipeline

**Expected output:** ~500 MB - 2 GB of data

---

## Features

- **Bulk download**: Efficient streaming API for large datasets
- **Comprehensive extraction**: Publications from both Swiss-Prot (curated) and TrEMBL (mapped)
- **Monthly granularity**: Time series data at month-level resolution
- **Multiple formats**: Long, wide, and summary formats for different analyses
- **Data validation**: Comprehensive quality checks and statistics
- **Error handling**: Automatic retry with exponential backoff
- **Progress tracking**: Real-time progress bars and logging

---

## Installation

### Requirements

- Python 3.8 or higher
- ~5 GB free disk space
- Internet connection

### Setup

1. **Clone/navigate to the uniprot directory:**

```bash
cd uniprot
```

2. **Create a virtual environment (recommended):**

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies:**

```bash
pip install -r requirements.txt
```

---

## Quick Start

### Run the full pipeline in sequence:

```bash
# Step 1: Download UniProt data (~30-60 minutes)
python scripts/01_download_proteins.py

# Step 2: Parse publications (~10-20 minutes)
python scripts/02_parse_publications.py

# Step 3: Aggregate by month (~5 minutes)
python scripts/03_aggregate_monthly.py

# Step 4: Validate data quality (~2 minutes)
python scripts/04_validate_data.py
```

### Run with custom directories:

```bash
# Download only reviewed (Swiss-Prot) proteins
python scripts/01_download_proteins.py --reviewed-only --output-dir data/raw

# Parse with custom directories
python scripts/02_parse_publications.py --input-dir data/raw --output-dir data/processed

# Aggregate with custom directories
python scripts/03_aggregate_monthly.py --input-dir data/processed --output-dir data/outputs

# Validate
python scripts/04_validate_data.py --processed-dir data/processed --outputs-dir data/outputs
```

---

## Pipeline Steps

### Step 1: Download Proteins (`01_download_proteins.py`)

Downloads human protein entries from UniProt REST API.

**Features:**
- Automatic retry with exponential backoff (2s, 4s, 8s, 16s)
- Streaming download for large files
- Downloads both reviewed (Swiss-Prot) and unreviewed (TrEMBL)
- Progress logging

**Output:**
- `data/raw/human_reviewed.json` - Swiss-Prot entries (~20 MB)
- `data/raw/human_unreviewed.json` - TrEMBL entries (~200+ MB)
- `data/raw/download_stats.json` - Download statistics

**Options:**
```bash
python scripts/01_download_proteins.py --help

Options:
  --reviewed-only     Only download reviewed (Swiss-Prot) entries
  --output-dir PATH   Output directory (default: data/raw)
```

**Expected output:**
```
Reviewed proteins: ~20,000 entries
Unreviewed proteins: ~150,000+ entries (if not using --reviewed-only)
```

---

### Step 2: Parse Publications (`02_parse_publications.py`)

Extracts publication data from downloaded JSON files.

**Features:**
- Extracts title, authors, journal, date, PMID, DOI
- Handles both curated and mapped citations
- Progress bars for large datasets
- Data cleaning and validation
- Multiple output formats (CSV + Parquet)

**Output:**
- `data/processed/all_publications.csv` - Full dataset in CSV
- `data/processed/all_publications.parquet` - Compressed Parquet format (recommended)
- `data/processed/parsing_statistics.json` - Parsing summary

**Options:**
```bash
python scripts/02_parse_publications.py --help

Options:
  --input-dir PATH    Input directory with downloaded files (default: data/raw)
  --output-dir PATH   Output directory (default: data/processed)
```

**Expected output:**
```
Total publications: ~50,000 - 200,000
Proteins with publications: ~15,000 - 20,000
Publications with PMID: ~90%+
Publications with month: ~60-80%
```

---

### Step 3: Aggregate Monthly (`03_aggregate_monthly.py`)

Creates monthly time series datasets from parsed publications.

**Features:**
- Monthly publication counts per protein
- Both long and wide format outputs
- Yearly aggregates
- Global time series (all proteins combined)
- Protein-level summary statistics

**Output:**
- `data/outputs/monthly_publications_long.csv` - Long format (protein-month rows)
- `data/outputs/monthly_publications_wide.csv` - Wide format (months as columns)
- `data/outputs/protein_summaries.csv` - Per-protein statistics
- `data/outputs/yearly_publications.csv` - Yearly aggregates
- `data/outputs/global_monthly_timeseries.csv` - Overall monthly trends
- `data/outputs/aggregation_statistics.json` - Aggregation summary

**Options:**
```bash
python scripts/03_aggregate_monthly.py --help

Options:
  --input-dir PATH    Input directory with parsed data (default: data/processed)
  --output-dir PATH   Output directory (default: data/outputs)
```

---

### Step 4: Validate Data (`04_validate_data.py`)

Performs comprehensive quality checks on all datasets.

**Features:**
- File existence checks
- Data completeness validation
- Year/month validity checks
- Distribution analysis
- Temporal consistency checks
- Sanity checks for well-known proteins

**Output:**
- `data/outputs/validation_results.json` - Detailed validation results
- Console output with check results
- Exit code 0 (pass) or 1 (fail)

**Options:**
```bash
python scripts/04_validate_data.py --help

Options:
  --processed-dir PATH  Processed data directory (default: data/processed)
  --outputs-dir PATH    Outputs directory (default: data/outputs)
  --output-file PATH    Validation results file (default: data/outputs/validation_results.json)
```

---

## Output Files

### Directory Structure

```
uniprot/
├── data/
│   ├── raw/                          # Downloaded UniProt JSON files
│   │   ├── human_reviewed.json       # Swiss-Prot entries
│   │   ├── human_unreviewed.json     # TrEMBL entries
│   │   └── download_stats.json       # Download metadata
│   │
│   ├── processed/                    # Parsed publication data
│   │   ├── all_publications.csv      # All publications (CSV)
│   │   ├── all_publications.parquet  # All publications (Parquet, recommended)
│   │   └── parsing_statistics.json   # Parsing summary
│   │
│   └── outputs/                      # Aggregated datasets
│       ├── monthly_publications_long.csv      # Monthly time series (long)
│       ├── monthly_publications_wide.csv      # Monthly time series (wide)
│       ├── protein_summaries.csv              # Protein-level statistics
│       ├── yearly_publications.csv            # Yearly aggregates
│       ├── global_monthly_timeseries.csv      # Global trends
│       ├── aggregation_statistics.json        # Aggregation summary
│       └── validation_results.json            # Quality check results
```

### Key Output Files

#### 1. `all_publications.csv` / `.parquet`

Full parsed publication dataset.

**Columns:**
- `accession`: UniProt ID (e.g., P04637)
- `gene_name`: Primary gene name (e.g., TP53)
- `protein_name`: Full protein name
- `reviewed`: Boolean - Swiss-Prot (True) or TrEMBL (False)
- `pmid`: PubMed ID
- `doi`: DOI identifier
- `title`: Publication title
- `publication_date`: Full date string (YYYY-MM or YYYY)
- `year`: Publication year (integer)
- `month`: Publication month (1-12, may be null)
- `journal`: Journal name
- `citation_type`: 'curated' or 'mapped'
- `categories`: Annotation categories
- `year_month`: Formatted year-month string

**Use cases:**
- Raw publication analysis
- PMID lookups
- Journal distribution studies
- Citation type comparisons

---

#### 2. `monthly_publications_long.csv`

Monthly publication counts in long format (one row per protein-month).

**Columns:**
- `accession`: UniProt ID
- `gene_name`: Gene name
- `protein_name`: Protein name
- `year_month`: Date (YYYY-MM-DD format, always first of month)
- `total_pubs`: Total publications that month
- `curated_pubs`: Curated publications
- `mapped_pubs`: Mapped publications

**Use cases:**
- Time series analysis
- Trend detection
- Protein activity tracking
- Easy filtering and grouping

---

#### 3. `monthly_publications_wide.csv`

Monthly publication counts in wide format (months as columns).

**Columns:**
- `accession`, `gene_name`, `protein_name`
- One column per month: `2020-01`, `2020-02`, etc.

**Use cases:**
- Matrix operations
- Correlation analysis
- Heatmaps
- Direct comparison across proteins

---

#### 4. `protein_summaries.csv`

Summary statistics for each protein.

**Columns:**
- `accession`, `gene_name`, `protein_name`
- `total_pubs`: Total publications
- `first_pub_year`: Year of first publication
- `last_pub_year`: Year of last publication
- `is_reviewed`: Swiss-Prot or TrEMBL
- `first_pub_date`: First publication date (with month)
- `last_pub_date`: Last publication date (with month)
- `pubs_with_month`: Publications with month information
- `avg_monthly_pubs`: Average publications per month
- `peak_month`: Month with most publications
- `peak_month_pubs`: Publication count in peak month

**Use cases:**
- Protein prioritization
- Research activity assessment
- Historical analysis
- Quick lookups

---

#### 5. `yearly_publications.csv`

Yearly publication counts per protein.

**Columns:**
- `accession`, `gene_name`, `protein_name`
- `year`: Publication year
- `total_pubs`: Total publications
- `curated_pubs`: Curated publications
- `mapped_pubs`: Mapped publications

**Use cases:**
- Annual trend analysis
- Longer-term patterns
- Year-over-year comparisons

---

#### 6. `global_monthly_timeseries.csv`

Global publication trends across all proteins by month.

**Columns:**
- `year_month`: Date
- `total_pubs`: Total publications across all proteins
- `curated_pubs`: Total curated publications
- `mapped_pubs`: Total mapped publications
- `num_proteins`: Number of proteins with publications that month

**Use cases:**
- Overall research activity trends
- Community-wide publication patterns
- Seasonal effects
- Data quality assessment

---

## Data Quality

### Expected Completeness

Based on UniProt data structure:

- **Publications with PMID:** 85-95%
- **Publications with month information:** 60-80%
- **Publications with DOI:** 40-60%
- **Reviewed (Swiss-Prot) proteins:** ~20,000 (100% have some annotation)
- **Unreviewed (TrEMBL) proteins:** ~150,000+ (variable annotation quality)

### Validation Checks

The validation script (`04_validate_data.py`) performs:

1. **File existence checks** - All expected files present
2. **Column completeness** - Required fields exist
3. **Data validity** - Years 1900-present, months 1-12
4. **Format validation** - UniProt accession format
5. **Distribution checks** - Reasonable publication counts
6. **Sanity checks** - Well-known proteins present (TP53, BRCA1, etc.)
7. **Temporal consistency** - Monthly totals match yearly totals
8. **No negative counts** - All counts ≥ 0

### Quality Metrics

To assess your extraction quality:

```bash
# Check validation results
cat data/outputs/validation_results.json

# View parsing statistics
cat data/processed/parsing_statistics.json

# View aggregation statistics
cat data/outputs/aggregation_statistics.json
```

---

## Known Limitations

### 1. Incomplete Month Information

**Issue:** ~20-40% of publications may only have year, not month

**Impact:** These publications are excluded from monthly time series but included in yearly aggregates

**Mitigation:** Use yearly aggregates for complete coverage; monthly for higher-resolution analysis

### 2. Update Lag

**Issue:** UniProt updates every 8 weeks

**Impact:** Most recent 2 months may have incomplete data

**Recommendation:** Note extraction date when analyzing recent trends

### 3. Citation Coverage

**Issue:** UniProt is conservative - not all papers mentioning a protein are included

**What's included:**
- Papers that functionally characterize proteins
- Papers that provide new sequences
- Papers that describe PTMs, interactions, etc.

**What's excluded:**
- Papers that only mention the protein in passing
- High-throughput screens without specific characterization

**Alternative:** For broader coverage, consider PubTator or Europe PMC

### 4. Historical Data

**Issue:** Older publications (pre-1980) may have incomplete metadata

**Impact:** Early date ranges may be less accurate

### 5. Non-human Citations

**Issue:** Papers may study orthologs or homologs from other species

**Impact:** Citation count includes cross-species research

---

## Troubleshooting

### Download Issues

**Problem:** Download fails or times out

**Solutions:**
1. Check internet connection
2. Script automatically retries with exponential backoff (4 attempts)
3. Try `--reviewed-only` flag for smaller download
4. Check UniProt API status: https://www.uniprot.org/help/api

**Problem:** Downloaded file is empty or corrupted

**Solutions:**
1. Delete the file and re-run download script
2. Check available disk space
3. Verify network stability during download

### Parsing Issues

**Problem:** "File not found" error

**Solution:** Ensure download completed successfully; check `data/raw/` directory

**Problem:** Parsing very slow

**Solution:** This is normal for large datasets; expect 10-20 minutes for full dataset

### Memory Issues

**Problem:** Out of memory errors

**Solutions:**
1. Close other applications
2. Process reviewed proteins only (`--reviewed-only`)
3. Increase system swap space
4. Use a machine with more RAM (8GB+ recommended)

### Validation Failures

**Problem:** Validation script reports failures

**Solutions:**
1. Review specific errors in validation output
2. Check if files were created by previous steps
3. Verify data integrity of downloaded files
4. Re-run pipeline from failing step

---

## Advanced Usage

### Processing Subsets

To test the pipeline on a smaller dataset:

```bash
# Download only reviewed proteins
python scripts/01_download_proteins.py --reviewed-only

# Then run remaining steps normally
python scripts/02_parse_publications.py
python scripts/03_aggregate_monthly.py
python scripts/04_validate_data.py
```

### Custom Analysis

After running the pipeline, use the output files for custom analysis:

```python
import pandas as pd

# Load monthly time series
monthly = pd.read_csv('data/outputs/monthly_publications_long.csv')
monthly['year_month'] = pd.to_datetime(monthly['year_month'])

# Analyze specific protein
tp53 = monthly[monthly['gene_name'] == 'TP53']
print(f"TP53 publications per month: {tp53['total_pubs'].mean():.2f}")

# Find proteins with increasing trends
# ... your analysis code ...
```

### Updating Data

To refresh data with latest UniProt release:

```bash
# Remove old raw data
rm -rf data/raw/*.json

# Re-run pipeline
python scripts/01_download_proteins.py
python scripts/02_parse_publications.py
python scripts/03_aggregate_monthly.py
python scripts/04_validate_data.py
```

---

## Citation

If you use this pipeline in your research, please cite:

**UniProt:**
> UniProt Consortium. "UniProt: the Universal Protein Knowledgebase in 2023."
> *Nucleic Acids Research* (2023). https://doi.org/10.1093/nar/gkac1052

**UniProt REST API:**
> Nightingale A, et al. "The UniProt REST API: Update and new features."
> *Nucleic Acids Research* (2024). https://doi.org/10.1093/nar/gkaf394

---

## Support

For issues specific to this pipeline:
- Check the troubleshooting section above
- Review log files: `download_proteins.log`, `parse_publications.log`, etc.
- Examine validation results: `data/outputs/validation_results.json`

For UniProt-specific questions:
- UniProt Help: https://www.uniprot.org/help
- UniProt API Documentation: https://www.uniprot.org/help/api

---

## License

This pipeline is provided as-is for research purposes. UniProt data is freely available under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

---

## Changelog

**Version 1.0 (October 2025)**
- Initial release
- Support for bulk download via streaming API
- Monthly and yearly aggregation
- Comprehensive validation suite
- Multiple output formats

---

**Happy analyzing!** If you discover interesting trends in protein research, we'd love to hear about them.
