# AlphaFold Impact Analysis: Complete Data Processing Pipeline

## Overview

This repository contains a comprehensive data processing pipeline to analyze AlphaFold's impact on biomedical research from 2020-2023. The pipeline transforms raw PubTator data (163M disease mentions, 72M gene mentions) into a balanced Gene-Month panel ready for regression analysis.

## Pipeline Architecture

```
Raw Data (12.5GB)
    ↓
[Raw Cleaning] → filter_pubtator_files.py
    ↓
Cleaned Data (4.5GB)
    ↓
[Phase 1] → 01_disease_filter.py
    ↓
Disease-Relevant PMIDs (21.8M)
    ↓
[Phase 2] → 02_gene_intersection.py
    ↓
Gene-Disease Intersections (33.6M)
    ↓
[Phase 3] → 03_temporal_filter.py
    ↓
AlphaFold Era Data (13.7M)
    ↓
[Phase 2A] → 02A_author_novelty.py (parallel)
    ↓
Author-Enriched Data (13.7M)
    ↓
[Phase 4] → 04_citation_enrichment.py
    ↓
Citation-Enriched Data (13.7M)
    ↓
[Phase 5] → 05_panel_construction.py
    ↓
Gene-Month Panel (945K obs)
    ↓
[Fix] → 07_fix_unique_diseases_simple.py
    ↓
Disease-Enriched Panel (945K obs)
    ↓
[Enhancement] → 08_add_gene_groups.py
    ↓
Final Analysis Dataset (945K obs)
```

## Script Flowchart

### Core Pipeline Scripts

1. **`scripts/raw_cleaning/filter_pubtator_files.py`**
   - **Input**: Original PubTator files (12.5GB)
   - **Output**: `cleaned/disease_pmid_mesh.txt`, `cleaned/gene_pmid_id.txt` (4.5GB)
   - **Purpose**: Extract only pmid+mesh_id and pmid+gene_id pairs
   - **Key Decision**: Preserve all data while reducing file size by 64%

2. **`scripts/data_processing/01_disease_filter.py`**
   - **Input**: `cleaned/disease_pmid_mesh.txt` (163M records)
   - **Output**: `processed/disease_relevant_pmids.parquet` (21.8M PMIDs)
   - **Purpose**: Filter to disease-relevant publications using MESH depth 3/4 + OMIM
   - **Key Decision**: Focus on clinical diseases, exclude overly broad categories

3. **`scripts/data_processing/02_gene_intersection.py`**
   - **Input**: `cleaned/gene_pmid_id.txt` + disease PMIDs + master protein dataset
   - **Output**: `processed/gene_disease_intersect.parquet` (33.6M records)
   - **Purpose**: Intersect gene mentions with disease publications and master gene list
   - **Key Decision**: Drop semicolon-separated gene IDs (310K records, 0.43%)

4. **`scripts/data_processing/03_temporal_filter.py`**
   - **Input**: Gene-disease intersections + publication dates
   - **Output**: `processed/gene_disease_temporal.parquet` (13.7M records)
   - **Purpose**: Filter to AlphaFold era (2020-2023) and add temporal variables
   - **Key Decision**: Focus on 2020-2023 as critical AlphaFold impact period

5. **`scripts/data_processing/02A_author_novelty.py`** (Parallel Track)
   - **Input**: Temporal data + OpenAlex author data
   - **Output**: `processed/author_novelty.parquet` (13.7M records)
   - **Purpose**: Add author information and newcomer/veteran flags
   - **Key Decision**: Use last author as proxy for research leadership

6. **`scripts/data_processing/04_citation_enrichment.py`**
   - **Input**: Author novelty data + iCite citation database
   - **Output**: `processed/citation_enriched.parquet` (13.7M records)
   - **Purpose**: Add citation counts and multi-temporal quality percentiles
   - **Key Decision**: Create 9 quality flags (annual/quarterly/bi-monthly × top 25%/10%/5%)

7. **`scripts/data_processing/05_panel_construction.py`**
   - **Input**: Citation-enriched data + protein metadata
   - **Output**: `final/final_balanced_panel.parquet` (945K obs)
   - **Purpose**: Aggregate to gene-month level and create balanced panel
   - **Key Decision**: Perfect balance (19,690 genes × 48 months = 945,120 observations)

8. **`scripts/data_processing/07_fix_unique_diseases_simple.py`** (Final Fix)
   - **Input**: Balanced panel + unique disease files
   - **Output**: `final/final_panel_COMPLETE.parquet/.dta` (945K obs)
   - **Purpose**: Merge unique disease innovation metrics
   - **Key Decision**: Simple gene_id + ym merge using proper time mapping

9. **`scripts/data_processing/08_add_gene_groups.py`** (Gene Groups Enhancement)
   - **Input**: Disease-enriched panel + genes_groups.dta
   - **Output**: `final/final_panel_COMPLETE.parquet` (945K obs with gene_group)
   - **Purpose**: Add functional gene group classifications
   - **Key Decision**: Left join to preserve panel structure, 73.7% coverage achieved

### Diagnostic/Fix Scripts (Not in Main Pipeline)

- `06_fix_final.py`: Early attempt at fixing time variables
- `06_fix_unique_diseases_and_stata_export.py`: Complex time mapping approach (abandoned)

## Key Data Processing Decisions

### 1. Disease Filtering Strategy
- **Included**: MESH depth 3/4 diseases + OMIM diseases
- **Excluded**: Overly broad categories (depth 1/2), overly specific subcategories (depth 5+)
- **Rationale**: Focus on clinically meaningful disease categories
- **Impact**: 60.98% retention (99.4M → 163M mentions)

### 2. Gene ID Handling
- **Challenge**: Semicolon-separated gene IDs (e.g., "9260;3783750") 
- **Decision**: Drop these entries entirely
- **Rationale**: Cannot reliably split and maintain data integrity
- **Impact**: Lost 310K records (0.43% of data)

### 3. Temporal Scope
- **Period**: 2020-2023 (AlphaFold announcement to maturation)
- **Rationale**: Captures critical impact period while maintaining manageable scope
- **Alternative**: Could extend to 2015+ for longer baseline (future work)

### 4. Author Novelty Definition
- **Focus**: Last author (research leader)
- **Baseline**: First appearance on specific gene within 2020-2023
- **Limitation**: May overestimate newcomer rate (no pre-2020 history)
- **Rationale**: Within-period dynamics still informative

### 5. Citation Quality Metrics
- **Multi-temporal**: Annual, quarterly, bi-monthly percentiles
- **Thresholds**: Top 25%, 10%, 5% within each time period
- **Rationale**: Enables robust temporal analysis across different granularities

### 6. Panel Construction
- **Structure**: Balanced panel (every gene × every month)
- **Missing Data**: Fill with zeros for inactive gene-months
- **Rationale**: Standard panel format for econometric analysis

### 7. Unique Disease Innovation
- **Challenge**: Different time numbering systems (our YYYYMM vs their sequential)
- **Solution**: Map via ym_date column in unique disease files
- **Final Approach**: Direct gene_id + ym merge after proper time mapping

### 8. Gene Group Classifications
- **Source**: Gene groups classification file (40,704 records)
- **Coverage**: 73.7% of panel genes have functional group assignments
- **Multiple Groups**: Some genes have multiple classifications (e.g., "442|454")
- **Decision**: Left join to preserve panel structure, enabling functional heterogeneity analysis

## Data Flow Summary

| Stage | Input Size | Output Size | Retention | Key Transformation |
|-------|------------|-------------|-----------|-------------------|
| Raw Cleaning | 12.5GB | 4.5GB | 64% | Extract pmid+id pairs only |
| Disease Filter | 163M mentions | 99M mentions | 61% | MESH depth 3/4 + OMIM filter |
| Gene Intersection | 72M genes | 34M genes | 47% | Disease relevance + master gene filter |
| Temporal Filter | 34M records | 14M records | 41% | 2020-2023 filter + temporal variables |
| Author Novelty | 14M records | 14M records | 100% | Add author metadata + novelty flags |
| Citation Enrichment | 14M records | 14M records | 100% | Add citations + quality percentiles |
| Panel Construction | 14M records | 945K obs | Aggregation | Gene-month aggregation + balance |
| Disease Innovation | 945K obs | 945K obs | 100% | Add unique disease metrics |
| Gene Groups | 945K obs | 945K obs | 100% | Add functional classifications (73.7% coverage) |

## Final Dataset Structure

**`final/final_panel_COMPLETE.parquet` (945,120 observations)**

### Panel Structure
- **Genes**: 19,690 (from master protein dataset)
- **Time Period**: 48 months (January 2020 - December 2023)
- **Balance**: Perfect (every gene × every month)
- **Activity Rate**: 75.3% of gene-months have publications

### Variables by Category

#### Identification & Time
- `gene_id`: Gene identifier (int64)
- `ym`: Year-month YYYYMM format (202001-202312)
- `year`, `month`: Temporal components
- `quarter`, `year_quarter`: Quarterly aggregation  
- `bimonth`, `year_bimonth`: Bi-monthly aggregation
- `ym_seq`: Sequential months 1-48 (for Stata)

#### Publication Metrics
- `n_papers`: Total papers per gene-month
- **Citation Quality (Annual)**: `n_top25_y`, `n_top10_y`, `n_top05_y`
- **Citation Quality (Quarterly)**: `n_top25_q`, `n_top10_q`, `n_top05_q`  
- **Citation Quality (Bi-monthly)**: `n_top25_b2`, `n_top10_b2`, `n_top05_b2`

#### Author Dynamics
- `n_newcomer_papers`: Papers with first-time authors on gene
- `n_veteran_papers`: Papers with repeat authors on gene

#### Disease Innovation
- `unique_mesh_count`: Total unique disease associations (82.1M total)
- `new_mesh_count`: New disease associations (6.2M total)

#### Gene Metadata  
- `protein_id`: Protein identifier
- `gene_name`: Gene symbol
- `protein_existence`: Evidence level
- `average_plddt`: AlphaFold confidence score
- `gene_group`: Functional gene group classification (73.7% coverage)

## Quality Validation

### Data Integrity Checks
- ✅ **Perfect Balance**: 19,690 × 48 = 945,120 observations
- ✅ **Author Consistency**: veteran + newcomer = total papers (0 mismatches)
- ✅ **Citation Coverage**: 93.4% of papers have citation data
- ✅ **Disease Innovation**: 75.3% coverage for unique diseases, 70.2% for new diseases
- ✅ **Gene Metadata**: 100% coverage

### Key Statistics
- **Total Publications**: 13.7M papers analyzed
- **Active Gene-Months**: 711,806 (75.3%)
- **High-Citation Papers**: 5.4M papers in top percentiles
- **Unique Authors**: 607K last authors identified
- **Disease Associations**: 88.4M total disease-gene links

## Usage Instructions

### For Python Analysis
```python
import pandas as pd
df = pd.read_parquet('final/final_panel_COMPLETE.parquet')

# Create AlphaFold treatment indicators
df['post_alphafold'] = (df['year'] >= 2021).astype(int)
df['post_announcement'] = (df['year'] >= 2020.5).astype(int)
```

### For Stata Analysis  
```stata
use "final/final_panel_COMPLETE.dta", clear
tsset gene_id ym_seq

* Create treatment indicators
gen post_alphafold = (year >= 2021)
gen post_announcement = (year >= 2020.5)
```

## Research Questions Enabled

1. **Research Activity**: Did AlphaFold increase publications per gene?
2. **Citation Quality**: Did AlphaFold affect high-impact research patterns?
3. **Author Entry**: How did AlphaFold impact newcomer researcher patterns?
4. **Disease Innovation**: Did AlphaFold accelerate unique disease associations?
5. **Temporal Dynamics**: Do effects vary across quarterly/bi-monthly granularity?
6. **Heterogeneity**: Are effects larger for high-confidence AlphaFold proteins?
7. **Functional Analysis**: Do AlphaFold effects vary by gene functional groups?

## Computational Requirements

### Processing Time
- Raw cleaning: ~5 minutes
- Full pipeline: ~30 minutes total
- Memory peak: ~8GB RAM
- Storage: ~15GB total (including intermediates)

### File Sizes
- Raw PubTator: 12.5GB
- Intermediate files: ~6GB  
- Final dataset: 7.7MB (Parquet), 180MB (Stata)

## Reproducibility Notes

- All scripts use `/usr/bin/python3` for consistency
- Processing statistics saved at each phase
- Parquet format ensures consistent data types
- Random seeds set where applicable
- Complete dependency documentation in requirements.txt

## Future Extensions

1. **Temporal Extension**: Include 2015-2019 for full author histories
2. **International Scope**: Include non-English publications
3. **Author Networks**: Full author team analysis (not just last author)
4. **Protein Families**: Group analysis by protein families/functions
5. **Citation Networks**: Include citing paper analysis
6. **Field Analysis**: Subject area classification of papers

---

*Dataset prepared for AlphaFold impact analysis. Ready for causal inference using difference-in-differences, event studies, and heterogeneous treatment effects models.*