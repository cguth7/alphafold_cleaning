# AlphaFold Impact Analysis: Results Summary

## 1. Matching Methodology

### 1.1 Treatment Definition
- **Treated**: Genes with `num_deposits == 0` (no prior PDB structure before AlphaFold)
- **Control**: Genes with `num_deposits > 0` (had prior experimental structure)
- **Treatment date**: July 2021 (AlphaFold structure release)

### 1.2 Propensity Score Matching (PSM)

**Step 1: Estimate propensity scores**
```
P(treated=1 | X) = logistic regression on:
  - pre_mean (avg papers/month pre-treatment)
  - pre_std, pre_trend
  - pre_newcomer, pre_veteran papers
  - pre_mesh, pre_new_mesh
  - average_plddt
  - protein_existence dummies
```

**Step 2: Nearest neighbor matching WITH REPLACEMENT**
- For each treated gene, find the control with closest propensity score
- Same control can match multiple treated genes

**Step 3: Weight controls by match frequency**
- Treated genes: weight = 1
- Control genes: weight = min(times_matched, 100)

### 1.3 Important Limitation

The matching procedure allows extreme control reuse:
- Gene 340096 matched **1,882 times** (capped to weight=100)
- Gene 144097 matched **263 times** (capped to weight=100)

This means these controls "represent" the counterfactual for many treated genes. Ideally, matching would cap reuse during the matching step (not just weights), but we lack the original matching code.

---

## 2. Sample Sizes

### 2.1 Overall Matched Sample
| Group | N Genes |
|-------|---------|
| Treated | 12,186 |
| Control | 3,817 |
| **Total** | **16,003** |

### 2.2 By Intensity Bin (pre-treatment papers/month)
| Bin | Treated | Control | T:C Ratio |
|-----|---------|---------|-----------|
| 0-1 | 5,311 | 352 | 15.1:1 |
| 1-3 | 3,812 | 1,368 | 2.8:1 |
| 3-5 | 1,256 | 676 | 1.9:1 |
| **5-10** | **1,017** | **691** | **1.5:1** |
| 10-20 | 457 | 399 | 1.1:1 |
| 20+ | 333 | 331 | 1.0:1 |

### 2.3 By PLDDT Quartile
| Quartile | Cutoff | Treated | Control | T:C Ratio |
|----------|--------|---------|---------|-----------|
| Q1 (low) | <65.4 | 3,275 | 729 | 4.5:1 |
| Q2 | 65.4-76.5 | 3,081 | 924 | 3.3:1 |
| Q3 | 76.5-85.7 | 2,962 | 1,038 | 2.9:1 |
| Q4 (high) | >85.7 | 2,880 | 1,126 | 2.6:1 |

---

## 3. DiD Results with Standard Errors

### 3.1 Aggregate DiD (different time aggregations)
| Frequency | β | SE | p-value | Sig |
|-----------|---|-----|---------|-----|
| Monthly | -0.226 | 0.071 | 0.0015 | ** |
| Quarterly | -0.649 | 0.219 | 0.0030 | ** |
| Semester | -1.298 | 0.453 | 0.0041 | ** |

### 3.2 DiD by Intensity Bin (Monthly)
| Bin | β | SE | p-value | Sig | N_T | N_C |
|-----|---|-----|---------|-----|-----|-----|
| 0-1 | -0.040 | 0.026 | 0.128 | | 5,311 | 352 |
| 1-3 | -0.085 | 0.025 | 0.0007 | *** | 3,812 | 1,368 |
| 3-5 | -0.027 | 0.061 | 0.663 | | 1,256 | 676 |
| **5-10** | **+0.146** | **0.069** | **0.033** | * | 1,017 | 691 |
| 10-20 | +0.136 | 0.181 | 0.453 | | 457 | 399 |
| 20+ | -1.484 | 1.724 | 0.389 | | 333 | 331 |

### 3.3 DiD by PLDDT Quartile (Monthly)
| Quartile | β | SE | p-value | Sig | N_T | N_C |
|----------|---|-----|---------|-----|-----|-----|
| Q1 (low) | -0.111 | 0.062 | 0.075 | + | 3,275 | 729 |
| Q2 | -0.151 | 0.082 | 0.065 | + | 3,081 | 924 |
| Q3 | -0.413 | 0.250 | 0.098 | + | 2,962 | 1,038 |
| Q4 (high) | -0.336 | 0.156 | 0.031 | * | 2,880 | 1,126 |

---

## 4. Monthly Event Studies (for reference)

![Monthly ES Cap 100](figures/event_studies_cap100.png)

---

## 5. Quarterly Event Studies

### 5.1 All Intensity Bins
![Quarterly ES All Bins](figures/quarterly_es_all_bins.png)

### 5.2 Aggregate (from all matched genes)

| Quarter | β | SE | 95% CI | Sig |
|---------|---|-----|--------|-----|
| -6 | +2.26 | 0.47 | [+1.34, +3.17] | *** |
| -5 | +1.90 | 0.42 | [+1.07, +2.72] | *** |
| -4 | +1.38 | 0.29 | [+0.82, +1.94] | *** |
| -3 | +1.11 | 0.21 | [+0.70, +1.52] | *** |
| -2 | -0.18 | 0.12 | [-0.41, +0.06] | |
| -1 | 0.00 | (ref) | - | |
| 0 | +0.23 | 0.12 | [-0.01, +0.46] | |
| 1 | -0.59 | 0.14 | [-0.87, -0.31] | *** |
| ... | | | | |
| 9 | +3.07 | 0.54 | [+2.02, +4.13] | *** |

**Pre-treatment avg: +1.08 | Post-treatment avg: +0.43**

**Interpretation**: Strong pre-trend violation. Treated genes were trending DOWN toward controls before treatment. Post-treatment shows initial dip then recovery.

### 5.3 Intensity Bin 5-10 (Best Identified)

| Quarter | β | SE | 95% CI | Sig |
|---------|---|-----|--------|-----|
| -6 | +0.45 | 0.36 | [-0.26, +1.16] | |
| -5 | +0.45 | 0.37 | [-0.28, +1.18] | |
| -4 | +0.31 | 0.39 | [-0.46, +1.09] | |
| -3 | +0.53 | 0.38 | [-0.23, +1.28] | |
| -2 | +0.84 | 0.40 | [+0.05, +1.63] | * |
| -1 | 0.00 | (ref) | - | |
| 0 | +0.15 | 0.39 | [-0.61, +0.92] | |
| 1 | +0.19 | 0.40 | [-0.60, +0.98] | |
| ... | | | | |
| 5 | +1.67 | 0.45 | [+0.78, +2.55] | *** |
| 6 | +1.53 | 0.45 | [+0.64, +2.42] | *** |
| 7 | +1.12 | 0.41 | [+0.30, +1.93] | ** |
| 8 | +1.54 | 0.39 | [+0.77, +2.31] | *** |
| 9 | +1.50 | 0.37 | [+0.78, +2.21] | *** |

**Pre-treatment avg: +0.43 | Post-treatment avg: +0.87**

**Interpretation**: Pre-trends are flatter (though not perfectly zero). Clear positive effect emerging ~5 quarters post-treatment, sustained through end of sample.

### 5.4 PLDDT Quartiles

![Quarterly ES PLDDT](figures/quarterly_es_plddt.png)

| Quartile | Pre-avg | Post-avg | Pattern |
|----------|---------|----------|---------|
| Q1 (low) | +0.38 | +0.11 | Pre-trend, weak post |
| Q2 | +1.08 | +0.65 | Pre-trend, moderate post |
| Q3 | +1.86 | +0.63 | Strong pre-trend |
| Q4 (high) | +1.36 | +0.36 | Pre-trend, U-shaped post |

**Interpretation**: All PLDDT quartiles show pre-trend violations. No clear pattern by PLDDT - effects don't vary systematically with prediction quality.

---

## 6. Semester Event Studies

### 6.1 Aggregate
| Semester | β | SE | 95% CI | Sig |
|----------|---|-----|--------|-----|
| -3 | +4.33 | 0.91 | [+2.55, +6.11] | *** |
| -2 | +2.67 | 0.49 | [+1.70, +3.63] | *** |
| -1 | 0.00 | (ref) | - | |
| 0 | -0.18 | 0.19 | [-0.56, +0.20] | |
| 1 | -0.23 | 0.25 | [-0.73, +0.27] | |
| 2 | -0.69 | 0.43 | [-1.52, +0.15] | |
| 3 | +1.41 | 0.32 | [+0.79, +2.03] | *** |
| 4 | +4.86 | 0.80 | [+3.29, +6.42] | *** |

### 6.2 Intensity 5-10 Bin
| Semester | β | SE | 95% CI | Sig |
|----------|---|-----|--------|-----|
| -3 | +0.07 | 0.57 | [-1.04, +1.17] | |
| -2 | +0.00 | 0.57 | [-1.11, +1.11] | |
| -1 | 0.00 | (ref) | - | |
| 0 | -0.50 | 0.56 | [-1.59, +0.59] | |
| 1 | +0.07 | 0.59 | [-1.10, +1.23] | |
| 2 | +0.96 | 0.67 | [-0.36, +2.28] | |
| 3 | +1.81 | 0.69 | [+0.45, +3.16] | ** |
| 4 | +2.20 | 0.61 | [+1.01, +3.39] | *** |

**Key**: 5-10 bin shows flat pre-trends at semester level and significant positive effects in semesters 3-4 (2023).

---

## 7. Key Findings

### 7.1 What We Can Say
1. **Aggregate effect is negative/null** depending on specification
2. **Strong heterogeneity by pre-treatment intensity**:
   - Low-activity genes (0-3 papers/mo): null or negative
   - **5-10 papers/mo: positive effect (~+0.15 papers/mo, p=0.03)**
   - High-activity (20+): unstable/unreliable
3. **No clear heterogeneity by PLDDT** - all quartiles show similar patterns
4. **Pre-trends are violated at aggregate level** - PSM did not fully solve this

### 7.2 Most Credible Result
The **5-10 papers/month intensity bin** is the most credible:
- Reasonable T:C ratio (1.5:1)
- Flattest pre-trends
- Significant positive effect that builds over time
- Robust to weight cap choice (10 vs 100)

### 7.3 Interpretation
> AlphaFold appears to have **accelerated research on moderately-studied proteins** (5-10 papers/month) but did **not stimulate new research on understudied proteins** (0-3 papers/month). The aggregate negative effect is driven by pre-trend violations and poor matching in low-activity bins.

---

## 8. Files

### Data
- `data/panel_matched_full.parquet` - Monthly panel with weights
- `data/panel_quarterly.parquet` - Quarterly aggregated
- `data/panel_semester.parquet` - Semester aggregated
- `data/matched_pairs.parquet` - PSM matched pairs

### Figures
- `figures/quarterly_es_all_bins.png` - Quarterly ES, all 6 intensity bins
- `figures/quarterly_es_intensity.png` - Quarterly ES, key intensity bins
- `figures/quarterly_es_plddt.png` - Quarterly ES by PLDDT quartile
- `figures/event_studies_cap100.png` - Monthly ES (cap=100)

---

*Analysis run: December 2024*
*Weight cap: 100*
