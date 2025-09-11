# AlphaFold Impact Analysis

**Measuring AlphaFold's Impact on Biomedical Research (2020-2023)**

This repository contains a complete data processing pipeline and analysis-ready dataset to study how AlphaFold affected biomedical research patterns, author entry, citation quality, and disease innovation.

## Quick Start

### Final Dataset
The analysis-ready dataset is available in `final/`:
- **`final_panel_COMPLETE.parquet`** (7.7 MB) - For Python analysis
- **Stata users**: Generate .dta file by running the processing scripts locally

```python
import pandas as pd
df = pd.read_parquet('final/final_panel_COMPLETE.parquet')
print(f"Dataset: {len(df):,} gene-month observations")
print(f"Genes: {df['gene_id'].nunique():,}, Months: {df['ym'].nunique()}")
```

### Key Research Questions
1. Did AlphaFold increase research activity per gene?
2. How did AlphaFold affect citation quality distributions?
3. Did AlphaFold impact newcomer author entry patterns?
4. Are effects heterogeneous across gene characteristics?
5. Do AlphaFold effects vary by functional gene groups?

## Dataset Overview

### Panel Structure
- **945,120 observations** = 19,690 genes × 48 months
- **Time Period**: January 2020 - December 2023 (AlphaFold era)
- **Perfect Balance**: Every gene has data for every month
- **Activity Rate**: 75.3% of gene-months have publications

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

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Contact

For questions about the dataset or methodology, please open an issue or contact [email].

---

**Ready for immediate analysis of AlphaFold's transformative impact on biomedical research!** 