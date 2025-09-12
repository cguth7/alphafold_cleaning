# AlphaFold Impact Analysis

**Measuring AlphaFold's Impact on Biomedical Research (2020-2023)**

This repository contains a complete data processing pipeline and analysis-ready dataset to study how AlphaFold affected biomedical research patterns, author entry, citation quality, and disease innovation.

## Quick Start

### Final Dataset
The analysis-ready datasets are available locally in `final/` (not tracked in git due to size):

**Recommended for Analysis:**
- **`final_panel_CLEAN.dta`** (172.6 MB) - Clean Stata format with essential columns only

**Full Datasets:**
- **`final_panel_COMPLETE_WITH_MISSING.dta`** (270.8 MB) - Complete with missing proteins flagged  
- **`final_panel_WITH_DEPOSITS.parquet`** (9.7 MB) - Python format with PDB deposits data
- **`final_panel_WITH_GROUPS.parquet`** (209.1 MB) - Complete with gene groups

**Archive:**
- **`final/old/`** - Previous versions moved for reference

```python
import pandas as pd
df = pd.read_stata('final/final_panel_CLEAN.dta')  # Recommended for analysis
print(f"Dataset: {len(df):,} gene-month observations")
print(f"Genes: {df['gene_id'].nunique():,}, Months: {df['ym'].nunique()}")
```

> **Note**: Data files are stored locally and not tracked in git due to GitHub's 100MB file limit. Run the processing scripts to generate the datasets.

### Key Research Questions
1. Did AlphaFold increase research activity per gene?
2. How did AlphaFold affect citation quality distributions?
3. Did AlphaFold impact newcomer author entry patterns?
4. Are effects heterogeneous across gene characteristics?
5. Do AlphaFold effects vary by functional gene groups?

## Dataset Overview

### Panel Structure
- **962,448 observations** = 20,051 genes × 48 months
- **Time Period**: January 2020 - December 2023 (AlphaFold era)
- **Perfect Balance**: Every gene has data for every month
- **Activity Rate**: 74.0% of gene-months have publications
- **Missing Flag**: 361 proteins (1.8%) flagged as `missing_from_literature = 1`

### Key Variables

#### Publication Metrics
- `n_papers`: Total papers per gene-month
- **Citation Quality**: Top 25%/10%/5% papers within year/quarter/bi-month
- **Author Dynamics**: Newcomer vs veteran author papers

#### Disease Innovation
- `unique_mesh_count`: Total unique disease associations (82M total)
- `new_mesh_count`: New disease associations (6.2M total)

#### Gene Information
- `gene_name`: Gene symbol
- `average_plddt`: AlphaFold confidence score
- `protein_existence`: Evidence level
- `gene_group`: Functional gene group classification (73.7% coverage)

## Data Processing Pipeline

The pipeline transforms 163M disease mentions and 72M gene mentions from PubTator into the final balanced panel:

```
Raw PubTator Data (12.5GB)
    ↓ [Raw Cleaning]
Cleaned Extracts (4.5GB)
    ↓ [Phase 1: Disease Filter]
Disease-Relevant PMIDs (21.8M)
    ↓ [Phase 2: Gene Intersection]  
Gene-Disease Links (33.6M)
    ↓ [Phase 3: Temporal Filter]
AlphaFold Era Data (13.7M)
    ↓ [Phase 2A: Author Novelty]
Author-Enriched (13.7M)
    ↓ [Phase 4: Citation Enrichment]
Citation-Enriched (13.7M) 
    ↓ [Phase 5: Panel Construction]
Gene-Month Panel (945K)
    ↓ [Final: Disease Innovation]
Disease-Enriched Panel (945K)
    ↓ [Gene Groups Addition]
Analysis-Ready Dataset (945K)
```

## Repository Structure

```
├── scripts/
│   ├── raw_cleaning/
│   │   └── filter_pubtator_files.py      # Extract pmid+id pairs
│   └── data_processing/
│       ├── 01_disease_filter.py          # Filter to disease-relevant PMIDs
│       ├── 02_gene_intersection.py       # Intersect genes with diseases
│       ├── 03_temporal_filter.py         # Filter to 2020-2023
│       ├── 02A_author_novelty.py         # Add author information
│       ├── 04_citation_enrichment.py     # Add citation metrics  
│       ├── 05_panel_construction.py      # Create balanced panel
│       ├── 07_fix_unique_diseases_simple.py  # Final disease merge
│       └── 08_add_gene_groups.py         # Add functional gene groups
├── processed/                            # Intermediate data files
├── final/                               # Analysis-ready datasets
├── PIPELINE_DOCUMENTATION.md           # Complete technical documentation
└── README.md                           # This file
```

## Key Processing Decisions

### Disease Scope
- **Included**: MESH depth 3/4 diseases + OMIM diseases
- **Rationale**: Focus on clinically meaningful categories
- **Impact**: 61% retention (163M → 99M mentions)

### Gene Coverage  
- **Source**: Master protein dataset (19,690 genes)
- **Coverage**: 99.2% of master genes found in literature
- **Focus**: Genes with structural information relevant to AlphaFold

### Temporal Scope
- **Period**: 2020-2023 (AlphaFold announcement to maturation)
- **Rationale**: Captures critical impact period
- **Events**: 2020 announcement, 2021 database launch, 2022-23 adoption

### Citation Quality
- **Metrics**: Multi-temporal percentiles (annual/quarterly/bi-monthly)
- **Thresholds**: Top 25%, 10%, 5% within time periods
- **Coverage**: 93.4% of papers have citation data

## Data Quality

### Validation Checks
- ✅ **Perfect Balance**: 19,690 × 48 = 945,120 observations
- ✅ **Author Consistency**: Veteran + newcomer = total papers (0 mismatches)
- ✅ **High Citation Coverage**: 5.4M high-quality papers identified
- ✅ **Disease Innovation**: 88.4M total disease-gene associations
- ✅ **Complete Metadata**: 100% gene information coverage

### Key Statistics
- **Total Publications**: 13.7M papers analyzed
- **Active Gene-Months**: 711,806 (75.3%)
- **Unique Authors**: 607K last authors
- **Newcomer Rate**: 84.4% (high research dynamism)

## Usage Examples

### Python Analysis
```python
import pandas as pd
import numpy as np

# Load data
df = pd.read_parquet('final/final_panel_COMPLETE.parquet')

# Create AlphaFold treatment indicators
df['post_alphafold_db'] = (df['year'] >= 2021).astype(int)
df['post_announcement'] = (df['year'] >= 2020.5).astype(int)

# Basic difference-in-differences setup
# Compare high vs low confidence proteins pre/post AlphaFold
df['high_confidence'] = (df['average_plddt'] > 70).astype(int)

# Example analysis: Effect on research activity
import statsmodels.api as sm
result = sm.ols(
    'n_papers ~ post_alphafold_db * high_confidence + gene_id + ym', 
    data=df
).fit()
```

### Stata Analysis

First, generate the .dta file from Python:
```python
import pandas as pd
df = pd.read_parquet('final/final_panel_COMPLETE.parquet')
df.to_stata('final/final_panel_COMPLETE.dta', write_index=False)
```

Then proceed with analysis:
```stata
use "final/final_panel_COMPLETE.dta", clear
tsset gene_id ym_seq

* Treatment indicators
gen post_alphafold = (year >= 2021)
gen high_confidence = (average_plddt > 70)

* Difference-in-differences
reghdfe n_papers i.post_alphafold##i.high_confidence, ///
    absorb(gene_id ym_seq) cluster(gene_id)
```

## Requirements

### Software
- Python 3.9+
- Required packages: `pandas`, `numpy`, `pyarrow`, `statsmodels`
- Stata 14+ (for .dta files - generate locally using the processing scripts)

### Hardware
- **Memory**: 8GB RAM recommended
- **Storage**: 15GB for full pipeline
- **Processing Time**: ~30 minutes for complete pipeline

## Installation

```bash
git clone https://github.com/cguth7/alphafold_cleaning.git
cd alphafold_cleaning
pip install -r requirements.txt
```

## Citations

If you use this dataset, please cite:

```bibtex
@dataset{alphafold_impact_2024,
  title={AlphaFold Impact Analysis Dataset: Gene-Month Panel 2020-2023},
  author={Guthmann, Max},
  year={2024},
  url={https://github.com/cguth7/alphafold_cleaning}
}
```

## Data Sources

- **PubTator**: Gene and disease mentions in biomedical literature
- **iCite**: Citation counts and metrics  
- **OpenAlex**: Author and institution information
- **AlphaFold DB**: Protein confidence scores

## Technical Notes

**Python:** All scripts use `/usr/bin/python3` for consistency across environments.

For complete technical documentation, see `PIPELINE_DOCUMENTATION.md`.

## Pipeline Enhancement Summary

This repository represents a major enhancement and correction of the original AlphaFold impact analysis pipeline. Key improvements and methodological corrections are documented below.

### 🎯 Major Enhancements Completed

#### 1. **Fixed Author Novelty Analysis (Section 2A)**
**Issue**: Original analysis only used 2020-2023 data to establish author baselines, incorrectly classifying many veteran authors as newcomers.

**Solution**: 
- Created extended temporal dataset (2015-2023) with 22.2M records
- Established proper author-gene baselines from 2015 onwards
- Corrected newcomer rate: 84.4% → 81.3% (more accurate veteran identification)

**Impact**: More accurate measurement of AlphaFold's effect on author entry patterns.

#### 2. **Integrated PDB Deposits Data**
**Enhancement**: Merged two deposit datasets into final panel:
- `deposits_dta.dta` (7,791 proteins) → `num_deposits` column  
- `deposits_monthly_dta.dta` (46,736 records) → `monthly_deposits` column

**Coverage**: 39.5% of proteins have deposit data, 1.5% of protein-months have monthly deposits.

#### 3. **Complete Protein Universe Recovery**
**Issue**: Original analysis missed 361 proteins from the master protein dataset.

**Solution**: Added missing proteins with dummy flag (`missing_from_literature = 1`).

**Result**: Complete balanced panel (20,051 proteins × 48 months = 962,448 records).

#### 4. **Enhanced Data Organization**
- **Clean Analysis File**: `final_panel_CLEAN.dta` with streamlined columns (24 vs 32)
- **Time Variable**: Single `ym` column (720-767 format) replacing 8 redundant time variables
- **Archive Structure**: Moved older versions to `final/old/` folder
- **Documentation**: Comprehensive column documentation and processing statistics

### ⚠️ Critical Methodological Differences

#### **MeSH Disease Levels: 3+4 vs 3+4+5**
**Key Finding**: The original Stata pipeline and this Python pipeline differ in disease scope:

- **Original Stata**: Uses MeSH depth levels **3, 4, AND 5** 
- **This Pipeline**: Uses MeSH depth levels **3 AND 4 only**

**Quantified Impact**:
- 88.6% of gene-months have identical `n_papers` counts
- 11.4% differ, with Stata consistently higher (+0.18 papers/gene-month average)  
- Total difference: ~121,000 additional paper counts in Stata version

**Implication**: Regression results will differ between pipelines due to this methodological choice. Including MeSH level 5 diseases captures more publications but may include less specific disease associations.

#### **Other Data Processing Decisions Affecting Results**

1. **Disease Filtering Strictness**
   - **Conservative** (this pipeline): Stricter OMIM + MeSH 3+4 only
   - **Inclusive** (Stata): OMIM + MeSH 3+4+5
   - **Impact**: ~8.5% difference in paper counts

2. **Gene Universe Definition**  
   - **Complete**: 20,051 proteins (includes 361 with no literature)
   - **Literature-Only**: 19,690 proteins (excludes proteins without publications)
   - **Impact**: Choice affects baseline and statistical power

3. **Author Baseline Period**
   - **Extended** (corrected): 2015-2023 for author history
   - **Limited** (original): 2020-2023 only
   - **Impact**: 3.1 percentage point reduction in newcomer rate

4. **Panel Balancing Strategy**
   - **Complete Fill**: All gene-months present (including zeros)
   - **Activity-Only**: Only gene-months with publications
   - **Impact**: Affects regression interpretation and standard errors

5. **Citation Percentile Calculations**
   - **Multi-Temporal**: Annual + Quarterly + Bi-monthly percentiles
   - **Single**: One time-based percentile only  
   - **Impact**: More granular citation quality measures

### 📊 Data Quality Validation

**Panel Balance**: ✅ Perfect (20,051 × 48 = 962,448)
**Time Coverage**: ✅ Complete (no missing months)
**Author Consistency**: ✅ Veteran + Newcomer = Total (zero mismatches)
**Citation Coverage**: ✅ 93.4% of papers have citation data
**Missing Data Handling**: ✅ Explicit flags for missing proteins

### 🔄 Reproducibility and Version Control

**Scripts**: All processing steps documented in numbered Python scripts
**Statistics**: Comprehensive processing statistics saved for each step
**Validation**: Built-in data quality checks and validation steps
**Paths**: Uses `/usr/bin/python3` for consistency across environments

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Contact

For questions about the dataset or methodology, please open an issue or contact [email].

---

**Ready for immediate analysis of AlphaFold's transformative impact on biomedical research!** 