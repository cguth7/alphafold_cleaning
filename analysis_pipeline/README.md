# AlphaFold2 Impact Analysis Pipeline

Python implementation replicating Danilo Messinese's CEM + Event Study analysis.

## Quick Start

```bash
python analysis_pipeline/run_pipeline.py
```

## What This Pipeline Does

1. **CEM Matching** (`01_cem_matching.py`)
   - Treatment: genes with NO prior PDB deposits (`num_deposits == 0`)
   - Control: genes WITH prior PDB deposits
   - Matching on pre-treatment (before July 2021) publication patterns:
     - `pre_mean_papers`: mean publications
     - `pre_mean_newcom`: mean newcomer papers
     - `pre_mean_veteran`: mean veteran papers
     - `pre_mean_top10`: mean top-10% papers
     - `pre_sd_papers`: std dev of papers
   - NOTE: pLDDT is NOT used for matching (per Matteo's instruction)

2. **Semester Aggregation** (`02_semester_aggregation.py`)
   - Aggregates monthly data to semester level
   - Creates ΔY (first difference) variables
   - Creates asinh transformations

3. **Event Study** (`03_event_study.py`)
   - Runs DiD event studies with three specifications:
     - LEVELS: raw outcome levels
     - DY: first differences (ΔY) - **preferred specification**
     - DASINH: asinh differences
   - Generates coefficient CSVs and figures

4. **DOI Integration** (`04_integrate_dois.py`)
   - Prepares DOI list for Shi & Evans metrics

## Key Outputs

```
analysis_pipeline/
├── data/
│   ├── matched_genes.parquet      # Gene-level CEM weights
│   ├── matched_panel_monthly.parquet
│   ├── matched_panel_semester.parquet
│   └── dois_for_shi_evans.csv     # For Danilo's S&E analysis
├── outputs/
│   └── eventstudy_*.csv           # Event study coefficients
└── figures/
    └── eventstudy_*.png           # Event study plots
```

## Key Results (DY Specification)

| Outcome | Pre-period Mean | Post-period Mean |
|---------|-----------------|------------------|
| n_papers | ~0.1 | ~2.2 |
| n_newcomer_papers | ~0.1 | ~1.4 |
| n_veteran_papers | ~0.0 | ~0.8 |
| n_top10_y | ~0.2 | ~0.6 |

Pre-period coefficients near zero indicate **flat pre-trends** (parallel trends assumption holds).

## Methodology Notes

This replicates Danilo's "CEM enriched" approach:
- Coarsened Exact Matching on pre-treatment publication patterns
- Event studies on ΔY (changes) rather than levels
- The combination of matching + differencing is what produces flat pre-trends

## Running Individual Steps

```bash
python analysis_pipeline/run_pipeline.py --step 1  # Only CEM matching
python analysis_pipeline/run_pipeline.py --step 3  # Only event studies
python analysis_pipeline/run_pipeline.py --skip-doi  # Skip DOI integration
```
