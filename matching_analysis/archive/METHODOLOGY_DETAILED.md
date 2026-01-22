# Precise Step-by-Step Methodology: PSM + DiD + Event Studies

## Overview

We want to estimate the causal effect of AlphaFold's protein structure release (July 2021) on research output (papers per gene per month).

---

## PHASE 1: DATA PREPARATION

### Step 1.1: Raw Panel Data

**Source file**: `final/final_panel_COMPLETE_WITH_MISSING.parquet`

**Structure**: Gene × Month panel
- 19,950 unique genes
- 48 months (Jan 2020 - Dec 2023)
- ~958,000 rows total

**Key variables**:
| Variable | Description |
|----------|-------------|
| `gene_id` | Unique gene identifier |
| `ym_seq` | Month sequence (1-48, where 1=Jan 2020) |
| `n_papers` | Count of papers mentioning this gene in this month |
| `num_deposits` | Number of PDB structures deposited for this gene BEFORE AlphaFold |

### Step 1.2: Define Treatment

**Treatment definition**:
```
treated = 1  if num_deposits == 0  (gene had NO prior experimental structure)
treated = 0  if num_deposits > 0   (gene HAD prior experimental structure)
```

**Intuition**: For genes with no prior structure, AlphaFold provides genuinely NEW structural information. For genes with existing structures, AlphaFold is somewhat redundant.

**Counts**:
- Treated: 12,179 genes
- Control: 7,771 genes

### Step 1.3: Define Post Period

**Treatment date**: July 2021 = `ym_seq = 19`

```
post = 1  if ym_seq >= 19  (July 2021 onwards)
post = 0  if ym_seq < 19   (before July 2021)
```

**Timeline**:
- Pre-period: 18 months (Jan 2020 - June 2021)
- Post-period: 30 months (July 2021 - Dec 2023)

---

## PHASE 2: PROPENSITY SCORE MATCHING (PSM)

### Step 2.1: Compute Pre-Treatment Features

For each gene, calculate characteristics using ONLY pre-treatment data (months 1-18):

| Feature | Definition |
|---------|------------|
| `pre_mean` | Average papers/month in pre-period |
| `pre_std` | Std deviation of papers in pre-period |
| `pre_trend` | Linear slope of papers over pre-period |
| `pre_newcomer` | Papers by new authors |
| `pre_veteran` | Papers by established authors |
| `pre_mesh` | Unique MeSH topic count |
| `pre_new_mesh` | New MeSH topics |
| `average_plddt` | AlphaFold prediction confidence score |
| `protein_existence_*` | Dummies for protein characterization level |

**Output file**: `data/gene_features_for_matching.parquet` (19,950 rows, one per gene)

### Step 2.2: Estimate Propensity Scores

**Model**: Logistic regression

```
P(treated = 1 | X) = logistic(β₀ + β₁·pre_mean + β₂·pre_std + β₃·pre_trend + ...)
```

**What this estimates**: The probability that a gene has no prior PDB structure, given its observable characteristics.

**Output**: `pscore` column added to features file
- Treated genes: mean pscore = 0.713, range [0.00, 1.00]
- Control genes: mean pscore = 0.450, range [0.00, 1.00]

### Step 2.3: Nearest Neighbor Matching

**Algorithm**:
```
FOR each treated gene T:
    1. Get T's propensity score: pscore_T
    2. Search all control genes
    3. Find control C that minimizes |pscore_T - pscore_C|
    4. Record the pair (T, C) and the distance
    5. C remains available for future matches (WITH REPLACEMENT)
```

**Key property - WITH REPLACEMENT**: The same control gene can be matched to multiple treated genes. This is necessary because we have 12,179 treated genes but only 7,771 controls with good pscore overlap.

**Output file**: `data/matched_pairs.parquet`

| Column | Description |
|--------|-------------|
| `treated_gene_id` | The treated gene |
| `control_gene_id` | Its matched control |
| `match_distance` | Absolute pscore difference |

**Statistics**:
- 12,179 pairs created (one per treated gene)
- 3,817 unique control genes used
- Average control is matched to 3.2 treated genes

### Step 2.4: Control Reuse Distribution

Because matching is with replacement, some controls are heavily reused:

| Times matched | # Controls |
|---------------|------------|
| Exactly 1 | 1,789 |
| 2-10 times | 1,941 |
| 11-100 times | 85 |
| >100 times | 2 |

**Extreme case**: Gene 340096 is matched **1,882 times**
- This gene has pscore = 0.997 (very high)
- It's matched to almost every high-pscore treated gene
- pre_mean = 0.4 papers/month (low activity)

---

## PHASE 3: CONSTRUCTING THE REGRESSION SAMPLE

### Step 3.1: Filter to Matched Genes

From the full panel (19,950 genes), keep only:
- All treated genes that appear in matched_pairs (12,179 genes)
- All control genes that appear in matched_pairs (3,817 genes)

**Result**: 15,996 genes × 48 months = 767,808 rows

### Step 3.2: Calculate Weights

**Purpose**: Account for control reuse in matching

**Formula**:
```
weight(gene g) =
    1                           if g is treated
    min(times_matched(g), 10)   if g is control
```

**Why cap at 10?** Without capping, gene 340096 (matched 1,882 times) would have weight 1,882 and dominate the entire estimate. Capping at 10 limits any single control's influence.

**Weight distribution**:
- Treated genes: all have weight = 1
- Control genes: mean weight = 3.2, max = 10 (after cap)

### Step 3.3: Create Regression Variables

```python
panel['treated'] = (num_deposits == 0).astype(int)      # 1 for treated, 0 for control
panel['post'] = (ym_seq >= 19).astype(int)              # 1 for July 2021+, 0 before
panel['treat_post'] = panel['treated'] * panel['post']  # DiD interaction term
```

The `treat_post` variable:
- = 1 ONLY for treated genes in post-period
- = 0 for controls (any time) or treated genes pre-treatment

---

## PHASE 4: DIFFERENCE-IN-DIFFERENCES REGRESSION

### Step 4.1: Model Specification

```
n_papers_it = αᵢ + γₜ + β·(treat_post)ᵢₜ + εᵢₜ
```

**Components**:

| Term | Meaning | What it absorbs |
|------|---------|-----------------|
| `αᵢ` | Gene fixed effect | All time-invariant gene differences (e.g., "BRCA1 always gets more papers") |
| `γₜ` | Time fixed effect | All common shocks (e.g., "COVID reduced papers in April 2020") |
| `β` | **TREATMENT EFFECT** | The causal parameter we want |
| `εᵢₜ` | Error term | Idiosyncratic variation |

### Step 4.2: What β Measures

β answers: **"How much did papers/month change for treated genes relative to control genes, after AlphaFold, beyond what we'd expect from gene-specific and time-specific factors?"**

Mathematically:
```
β = [E(Y|treated,post) - E(Y|treated,pre)] - [E(Y|control,post) - E(Y|control,pre)]
```

This is the "difference in differences":
- First difference: treated before vs after
- Second difference: subtract control before vs after
- Removes common trends, isolates treatment effect

### Step 4.3: Estimation Details

**Estimator**: Weighted Least Squares (via PanelOLS)

**Weighting**: Each observation weighted by `weight` from Step 3.2

**Objective function**:
```
minimize Σᵢₜ weightᵢ × (n_papersᵢₜ - αᵢ - γₜ - β·treat_postᵢₜ)²
```

**Standard errors**: Clustered at gene level
- Allows for arbitrary correlation within a gene over time
- Accounts for serial correlation in errors

### Step 4.4: Aggregate Results

| Specification | β | SE | p-value |
|---------------|---|-----|---------|
| PSM weighted, capped at 10 | **-0.276** | 0.075 | 0.0002 |

**Interpretation**: Treated genes saw 0.28 fewer papers/month relative to controls after AlphaFold, controlling for gene and time effects.

---

## PHASE 5: HETEROGENEOUS EFFECTS BY INTENSITY BIN

### Step 5.1: Define Bins

Partition genes by pre-treatment research activity (`pre_mean`):

| Bin | Range (papers/month) |
|-----|---------------------|
| 0-1 | [0, 1) |
| 1-3 | [1, 3) |
| 3-5 | [3, 5) |
| 5-10 | [5, 10) |
| 10-20 | [10, 20) |
| 20+ | [20, ∞) |

### Step 5.2: Sample Sizes per Bin

| Bin | Treated | Control | T:C Ratio |
|-----|---------|---------|-----------|
| 0-1 | 5,311 | 352 | 15:1 |
| 1-3 | 3,812 | 1,368 | 2.8:1 |
| 3-5 | 1,256 | 676 | 1.9:1 |
| 5-10 | 1,017 | 691 | 1.5:1 |
| 10-20 | 457 | 399 | 1.1:1 |
| 20+ | 333 | 331 | 1.0:1 |

**Problem**: The 0-1 bin has 15:1 ratio - very few controls to match against.

### Step 5.3: Run Separate Regressions

**CRITICAL**: For each bin, we run a **SEPARATE** regression:

```python
for bin in ['0-1', '1-3', '3-5', '5-10', '10-20', '20+']:
    # Subset to genes in this bin
    bin_data = matched_panel[matched_panel['intensity_bin'] == bin]

    # Run the SAME DiD specification
    model = PanelOLS(
        bin_data['n_papers'],
        bin_data[['treat_post']],
        entity_effects=True,    # gene FE
        time_effects=True,      # time FE
        weights=bin_data['weight']
    )
    results[bin] = model.fit(cov_type='clustered', cluster_entity=True)
```

**What this means**: The 5-10 bin regression ONLY uses genes with pre_mean in [5, 10). It compares treated vs control genes WITHIN that activity level.

### Step 5.4: Results by Bin

| Bin | β | SE | p-value | Sig at 5%? |
|-----|---|-----|---------|------------|
| 0-1 | -0.060 | 0.018 | 0.001 | Yes (negative) |
| 1-3 | -0.089 | 0.021 | <0.001 | Yes (negative) |
| 3-5 | -0.031 | 0.062 | 0.62 | No |
| 5-10 | **+0.148** | 0.069 | **0.031** | **Yes (positive)** |
| 10-20 | +0.136 | 0.181 | 0.45 | No |
| 20+ | -1.484 | 1.724 | 0.39 | No |

---

## PHASE 6: EVENT STUDY

### Step 6.1: Purpose

The event study:
1. **Tests parallel trends**: Pre-treatment coefficients should be ~0
2. **Shows timing**: When does the effect appear?
3. **Visualizes dynamics**: Does effect grow/fade over time?

### Step 6.2: Create Relative Time Variable

```python
treatment_month = 19  # July 2021
panel['rel_month'] = panel['ym_seq'] - treatment_month
```

**Result**:
- rel_month = -18: January 2020
- rel_month = -1: June 2021 (reference period)
- rel_month = 0: July 2021 (treatment)
- rel_month = +29: December 2023

### Step 6.3: Create Event Study Dummies

For each relative month k (except reference month -1):

```python
panel[f'treat_m{k}'] = panel['treated'] * (panel['rel_month'] == k)
```

This creates 47 dummy variables:
- `treat_m-18`: = 1 if treated AND month is Jan 2020
- `treat_m-17`: = 1 if treated AND month is Feb 2020
- ... (skip -1 as reference)
- `treat_m0`: = 1 if treated AND month is July 2021
- ...
- `treat_m29`: = 1 if treated AND month is Dec 2023

### Step 6.4: Event Study Regression

```
n_papers_it = αᵢ + γₜ + Σₖ βₖ·(treat_mₖ)ᵢₜ + εᵢₜ
                       (k ≠ -1)
```

**What each βₖ measures**: The difference between treated and control at relative month k, relative to the difference at month -1 (the reference).

### Step 6.5: Interpreting Event Study Coefficients

**Pre-treatment (k < 0)**:
- If parallel trends holds: all βₖ ≈ 0
- If βₖ ≠ 0 before treatment: parallel trends violated, DiD is biased

**Post-treatment (k ≥ 0)**:
- βₖ > 0: treated genes doing better than controls (positive effect)
- βₖ < 0: treated genes doing worse than controls (negative effect)
- Pattern shows effect dynamics (immediate jump? gradual build?)

### Step 6.6: Event Study by Bin

Same logic as DiD: run SEPARATE event studies for each intensity bin.

```python
for bin in bins:
    bin_data = matched_panel[matched_panel['intensity_bin'] == bin]
    # Run event study regression on bin_data only
```

---

## PHASE 7: KEY RESULTS

### Aggregate Event Study

| Period | Avg coefficient | Interpretation |
|--------|-----------------|----------------|
| Pre-treatment (months -18 to -2) | **+0.37** | NOT ZERO - parallel trends violated |
| Post-treatment (months 0 to 29) | +0.07 | Near zero |

**Problem**: Pre-treatment coefficients are significantly positive and trending downward. This means PSM failed to achieve parallel trends at the aggregate level.

### Event Study by Bin

| Bin | Pre avg | Post avg | Parallel trends? | Effect? |
|-----|---------|----------|------------------|---------|
| 0-1 | -0.16 | -0.21 | Questionable | Negative |
| 1-3 | -0.23 | -0.31 | No (drifting) | Negative |
| 3-5 | -0.51 | -0.51 | Flat but negative | Null |
| **5-10** | **+0.09** | **+0.23** | **≈ Yes** | **Positive** |
| 10-20 | +0.18 | +0.31 | Noisy | Unclear |
| 20+ | +2.17 | +0.57 | No (unstable) | Unreliable |

### Best Identified Result: 5-10 Bin

- Pre-treatment coefficients fluctuate around zero (parallel trends approximately hold)
- Post-treatment coefficients gradually become positive
- Effect builds over time (months 16+ show significant positive coefficients)
- β = +0.15 papers/month, p = 0.03

---

## SUMMARY FLOWCHART

```
RAW DATA (gene × month panel)
         ↓
CALCULATE PRE-TREATMENT FEATURES (for each gene)
         ↓
ESTIMATE PROPENSITY SCORES (logistic regression)
         ↓
NEAREST NEIGHBOR MATCHING (with replacement)
         ↓
CREATE MATCHED SAMPLE + WEIGHTS
         ↓
    ┌────┴────┐
    ↓         ↓
   DiD      Event Study
    ↓         ↓
 Single β   β for each month
    ↓         ↓
    └────┬────┘
         ↓
REPEAT FOR EACH INTENSITY BIN (separate regressions)
```

---

# QUIZ: Test Your Understanding

Answer these questions to check if you understand the methodology:

## Question 1: Treatment Definition
What makes a gene "treated" in this analysis? Why is this a sensible treatment definition?

## Question 2: Propensity Score
A gene has pscore = 0.85. What does this number mean? Is this gene more likely to be treated or control?

## Question 3: Matching With Replacement
Control gene X is matched to 50 different treated genes.
- a) What weight does gene X get in the regression?
- b) Why do we cap weights at 10?
- c) If we didn't cap, what problem would arise?

## Question 4: The DiD Interaction
In the regression `n_papers = α_i + γ_t + β(treat_post) + ε`:
- a) When exactly does treat_post = 1?
- b) When does treat_post = 0?
- c) What does β estimate?

## Question 5: Fixed Effects
- a) What does the gene fixed effect αᵢ absorb?
- b) What does the time fixed effect γₜ absorb?
- c) Why are both necessary?

## Question 6: Event Study Pre-Trends
The aggregate event study shows pre-treatment coefficients averaging +0.37 (significantly different from zero).
- a) What does this tell us about parallel trends?
- b) What does this imply about the aggregate DiD estimate?

## Question 7: Bin-Specific Analysis
For the 5-10 papers/month bin, we get β = +0.15, p = 0.03.
- a) What sample is used for this regression? (All genes? Only some?)
- b) Why might this bin be more credible than the 0-1 bin?

## Question 8: Weight Interpretation
In weighted OLS, what does it mean mathematically for a control gene to have weight = 5?

## Question 9: Critical Thinking
The 0-1 bin has 5,311 treated genes but only 352 controls (15:1 ratio).
- a) Why is this problematic for the matching?
- b) What happens to the extreme control (gene 340096, matched 1882 times)?

## Question 10: Bottom Line
Based on this analysis, can we say AlphaFold caused more research? What's the most defensible claim?

---

*Answers provided upon request*
