# Processed Data Directory

This directory contains data files from each processing phase, optimized for analysis and saved in Parquet format.

## Phase 1: Disease-Relevant Publication Filter

### disease_relevant_pmids.parquet
**Source:** `cleaned/disease_pmid_mesh.txt` filtered by MESH depth 3/4 or OMIM diseases  
**Size:** 95.2 MB  
**Rows:** 21,844,614 unique PMIDs  
**Processing Script:** `scripts/data_processing/01_disease_filter.py`

**Columns:**
- `pmid` (int32): PubMed ID of disease-relevant publications

**Processing Summary:**
- **Input:** 163,082,572 disease mentions from PubTator
- **Filter:** MESH depth 3/4 diseases OR OMIM diseases  
- **Retention:** 99,440,126 mentions (60.98% kept)
- **Unique PMIDs:** 21,844,614 publications
- **Average mentions per paper:** 4.55 diseases per publication

**Quality Metrics:**
- OMIM diseases kept: 1,077,433 mentions
- MESH depth 3/4 kept: 98,362,693 mentions  
- Total file compression: 3.4GB → 95.2MB (97.2% reduction)

**Sample Excluded MESH IDs:**
- MESH:D013167 (ankylosing spondylitis - likely depth 5)
- MESH:D007669 (kidney diseases - likely depth 5)  
- MESH:D012559 (liver diseases - likely depth 5)

See `phase1_processing_stats.txt` for complete processing statistics.

## Phase 2: Gene-Disease Intersection

### gene_disease_intersect.parquet
**Source:** Gene mentions from `cleaned/gene_pmid_id.txt` filtered by disease-relevant PMIDs and master gene list  
**Size:** 118.0 MB  
**Rows:** 33,618,023 gene mentions  
**Processing Script:** `scripts/data_processing/02_gene_intersection.py`

**Columns:**
- `pmid` (int32): PubMed ID of disease-relevant publications
- `gene_id` (int64): Gene identifier from master protein dataset

**Processing Summary:**
- **Step 1: Disease Intersection**
  - Input: 72.6M gene mentions (after dropping 310K semicolon-separated IDs)
  - Disease-relevant: 59.6M mentions (82.2% retained)
  - Unique PMIDs: 6.7M → 5.6M publications
  - Gene retention: 49.2% (many genes only in non-disease papers)

- **Step 2: Master Gene Filter**
  - Master genes available: 19,950 from protein dataset
  - Genes matched: 19,783 (99.2% of master genes found)
  - Final mentions: 33.6M (56.4% of disease-relevant mentions)
  - **Key insight:** Only 0.6% of PubTator genes are in master protein dataset

**Quality Metrics:**
- Average mentions per paper: 5.99 genes per publication
- Data compression: 1.4GB → 118MB (92% reduction)  
- Semicolon-separated gene IDs dropped: 310,777 (0.43%)

**Important Notes:**
- ⚠️ **Low gene retention (0.6%)**: Most genes in PubTator are NOT in the master protein dataset
- This is expected - master dataset focuses on specific proteins with structural data
- Sample excluded genes: 64030, 24887, 24224 (likely valid genes not in AlphaFold scope)

See `phase2_processing_stats.txt` for detailed merge statistics.

## Phase 3: Temporal Filtering (AlphaFold Era)

### gene_disease_temporal.parquet
**Source:** Gene-disease intersections from Phase 2 with publication dates, filtered to 2020-2023  
**Size:** 57.2 MB  
**Rows:** 13,721,571 gene mentions  
**Processing Script:** `scripts/data_processing/03_temporal_filter.py`

**Columns:**
- `pmid` (int32): PubMed ID of disease-relevant publications
- `gene_id` (int64): Gene identifier from master protein dataset
- `year` (int16): Publication year  
- `month` (int8): Publication month
- `ym` (int64): Year-month combined (YYYYMM format)
- `quarter` (int8): Quarter within year (1-4)
- `year_quarter` (int64): Year-quarter combined (YYYYQ format)
- `bimonth` (int8): Bi-month within year (1-6)
- `year_bimonth` (int64): Year-bimonth combined (YYYYB format)

**Processing Summary:**
- **Step 1: Date Merging**
  - Gene-disease records with dates: 30.6M (90.9% merge success)
  - PMID merge rate: 80.2% (4.5M PMIDs matched)
  - Missing dates dropped: 16,243 PMIDs (0.24%)

- **Step 2: AlphaFold Era Filtering (2020-2023)**
  - Temporal retention: 44.9% of dated records  
  - Final publications: 1.45M PMIDs
  - Gene retention: 99.6% (19,690/19,779 genes)
  - Time focus: Critical AlphaFold impact period

**Temporal Distribution:**
- **2020:** 2.96M mentions (21.6%) - AlphaFold announcement
- **2021:** 3.80M mentions (27.7%) - AlphaFold DB launch  
- **2022:** 3.91M mentions (28.5%) - Peak adoption
- **2023:** 3.06M mentions (22.3%) - Maturation

**Quality Metrics:**
- Comprehensive time variables for quarterly/bi-monthly analysis
- Data compression: 384MB → 57MB (85% reduction)
- Perfect gene coverage retention across time filtering
- Ready for citation enrichment and panel construction

**Key Insights:**
- Strong publication growth 2020-2022 aligns with AlphaFold milestones
- Excellent gene retention shows master proteins well-represented in recent literature
- Date merge success indicates high-quality PMID coverage

See `phase3_processing_stats.txt` for detailed temporal statistics.

## Phase 2A: Author Novelty Analysis (Parallel Track)

### author_novelty.parquet
**Source:** Temporal gene-disease data merged with deduplicated author information  
**Size:** 132.5 MB  
**Rows:** 13,673,571 author-gene records  
**Processing Script:** `scripts/data_processing/02A_author_novelty.py`

**Columns:**
- All Phase 3 columns (pmid, gene_id, year, month, temporal variables)
- `last_author_id` (string): OpenAlex author identifier  
- `last_author_name` (string): Author name
- `last_institution_name` (string): Author's institution
- `last_institution_id` (string): OpenAlex institution identifier
- `newcomer_author` (int): 1 for first publication by author on gene, 0 otherwise
- `not_newcomer_author` (int): 1 for repeat publication by author on gene, 0 otherwise

**Processing Summary:**
- **Step 1: Author Merging**
  - Author merge success: 99.65% (13.7M → 13.7M records)
  - PMID coverage: 99.58% (1.45M → 1.44M PMIDs with authors)
  - Authors identified: 607,477 unique last authors

- **Step 2: Novelty Analysis**  
  - Author-gene pairs: 11.54M unique combinations
  - **Newcomer rate: 84.4%** (first-time author-gene publications)
  - **Veteran rate: 15.6%** (repeat author-gene publications)

**Key Insights:**
- **High newcomer rate (84.4%)** suggests dynamic research landscape during AlphaFold era
- Limited by 2020-2023 baseline (ideal would include 2015+ for complete author histories)
- Excellent author coverage (99.6%) shows comprehensive OpenAlex integration
- Ready for gene-month aggregation to track author entry patterns

**Usage in Analysis:**
- Aggregate newcomer/veteran counts per gene-month
- Track how AlphaFold affected new researcher entry into gene research
- Compare author novelty patterns across genes and time periods

**Important Note:**
- Baseline limited to 2020-2023 data (vs ideal 2015+)
- May overestimate newcomer rate for authors active pre-2020
- Still valuable for within-period author dynamics analysis

See `phase2A_processing_stats.txt` for detailed author merge statistics.

## Phase 4: Citation Enrichment

### citation_enriched.parquet
**Source:** Author novelty data enriched with iCite citation metrics and quality flags  
**Size:** 153.1 MB  
**Rows:** 13,673,571 citation-enriched records  
**Processing Script:** `scripts/data_processing/04_citation_enrichment.py`

**New Columns Added:**
- `icite_cites` (int): Number of citations from iCite database
- **Annual flags**: `hi25_year`, `hi10_year`, `hi05_year` (top 25%/10%/5% within year)
- **Quarterly flags**: `hi25_quarter`, `hi10_quarter`, `hi05_quarter` (top 25%/10%/5% within quarter)
- **Bi-monthly flags**: `hi25_bimonth`, `hi10_bimonth`, `hi05_bimonth` (top 25%/10%/5% within bi-month)

**Processing Summary:**
- **Citation Coverage**: 93.4% of PMIDs have citation data (1.35M/1.44M papers)
- **Citation Statistics**:
  - Mean citations: 24.8 per paper
  - Median citations: 11 per paper  
  - Range: 1-8,090 citations
  - Missing citations filled with 0

**High-Citation Paper Distribution:**
- **Annual**: 376K top-25% papers (26.1%), 147K top-10% (10.2%), 73K top-5% (5.1%)
- **Quarterly**: 381K top-25% papers (26.4%), 149K top-10% (10.3%), 74K top-5% (5.1%)  
- **Bi-monthly**: 379K top-25% papers (26.2%), 148K top-10% (10.2%), 74K top-5% (5.1%)

**Author Dynamics in High-Quality Research:**
- **92.3% of papers** have newcomer authors (first-time on gene)
- **35.1% of papers** have veteran authors (repeat on gene)
- Shows dynamic research landscape with many new entrants

**Key Features:**
- **Multi-granularity percentiles** enable fine-grained temporal analysis
- **Paper-level aggregation** ensures consistent citation flags across gene mentions
- **Comprehensive coverage** with excellent citation database integration
- **Ready for aggregation** at gene-month level with quality stratification

**Usage in Analysis:**
- Compare high-citation publication patterns before/after AlphaFold
- Analyze whether AlphaFold affected citation quality distributions
- Track newcomer vs veteran success in high-impact research
- Multi-temporal granularity supports robust temporal analysis

See `phase4_processing_stats.txt` for detailed citation merge statistics.

---

## Data Quality Notes

- All PMIDs are validated integers
- No missing values in output
- Parquet format provides efficient compression and fast loading
- Processing statistics preserved for reproducibility