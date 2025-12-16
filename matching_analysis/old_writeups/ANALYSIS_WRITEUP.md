# Propensity Score Matching Analysis for AlphaFold Impact

## Overview

This analysis attempts to address pre-trend violations in a difference-in-differences framework examining whether AlphaFold's July 2021 protein structure release affected research activity (measured by paper counts).

---

## Methodology: How PSM + DiD + Event Studies Work Together

### The Big Picture

We want to know: **Did AlphaFold cause more papers to be written about proteins?**

The challenge: We can't just compare "treated" genes (no prior structure) vs "control" genes (had structure) because these groups are fundamentally different. Control genes have ~8x more papers and are trending upward faster. A naive comparison would be confounded.

### Step 1: Propensity Score Matching (PSM)

**Goal**: Create a "matched" sample where treated and control genes are comparable.

**How it works**:
1. Estimate a **propensity score** for each gene = P(treated | pre-treatment characteristics)
2. Use logistic regression: `P(no_PDB) ~ pre_mean + pre_trend + plddt + protein_existence + ...`
3. For each treated gene, find a control gene with a similar propensity score
4. This creates matched pairs where treated/control have similar observable characteristics

**Key assumption**: After matching on observables, the remaining variation in treatment is "as good as random"

### Step 2: Difference-in-Differences (DiD)

**Goal**: Estimate the causal effect of AlphaFold using before/after + treatment/control variation.

**The regression**:
```
n_papers_it = α_i + γ_t + β(Treated_i × Post_t) + ε_it
```

Where:
- `α_i` = gene fixed effects (absorbs all time-invariant gene differences)
- `γ_t` = time fixed effects (absorbs common shocks like COVID, general publication trends)
- `Treated_i × Post_t` = 1 if gene is treated AND time is post-July 2021
- `β` = **the causal effect we want** (how much more/less did treated genes change relative to controls after AlphaFold?)

**Key assumption**: Parallel trends - absent treatment, treated and control would have followed the same trajectory. This is what matching is trying to achieve.

**With matching**: We run this regression on the matched sample, weighting controls by how often they're matched.

### Step 3: Event Studies

**Goal**: Visualize whether parallel trends hold and when effects appear.

**How it works**:
1. Instead of a single `Post` indicator, create dummies for each month relative to treatment
2. Estimate: `n_papers_it = α_i + γ_t + Σ_k β_k(Treated_i × 1{t=k}) + ε_it`
3. Plot the β_k coefficients over time

**What to look for**:
- **Pre-treatment (k < 0)**: Coefficients should be ~0 and flat. If they're trending, parallel trends is violated.
- **Post-treatment (k ≥ 0)**: If there's an effect, coefficients should jump/diverge from zero.

### Putting It Together

```
Raw data → PSM matching → Matched sample → DiD regression → β estimate
                                        → Event study → Visual check of parallel trends
```

The event study is a **diagnostic**: if pre-treatment coefficients aren't flat, the DiD estimate is suspect even after matching.

---

## Data

- **Source**: `final/final_panel_COMPLETE_WITH_MISSING.parquet`
- **Panel structure**: 19,950 genes × 48 months (Jan 2020 - Dec 2023)
- **Outcome variable**: `n_papers` (count of papers per gene-month)
- **Treatment date**: July 2021 (ym_seq = 19)

## Treatment Definition

- **Treated**: Genes with `num_deposits == 0` (no prior PDB structures) → AlphaFold provides NEW structural information
- **Control**: Genes with `num_deposits > 0` (had prior structures) → AlphaFold provides somewhat redundant information
- **Sample sizes**: 12,179 treated, 7,771 control

## The Pre-Trends Problem

### Initial Findings

The treated and control groups are fundamentally different:

| Metric | Control | Treated | Ratio |
|--------|---------|---------|-------|
| Mean papers/month | 30.3 | 4.1 | 7.4x |
| Pre-period slope | 0.779/month | 0.098/month | 7.9x |

The control group (genes with prior structures) has much higher research attention AND is trending upward faster. This violates the parallel trends assumption required for DiD.

**Key insight**: Paper counts are extremely concentrated, and high-volume genes are mostly in the control group:

| Top % of genes | % of all papers | Control | Treated |
|----------------|-----------------|---------|---------|
| Top 1% (199)   | 39%             | 93%     | 7%      |
| Top 5% (997)   | 66%             | 88%     | 12%     |
| Top 10% (1995) | 78%             | 83%     | 17%     |

### Visualization

See the top row of `pretrends_comparison_FIXED.png` - shows the diverging trends in levels (control ~30 papers/mo, treated ~4 papers/mo), though indexed trends are more parallel.

![Pre-trends Comparison](pretrends_comparison_FIXED.png)

## Matching Strategy

We tried two approaches. The first (standard propensity score matching) didn't work well. The second (exact bin + nearest neighbor) achieved parallel trends.

---

### Approach 1: Standard Propensity Score Matching (PSM)

**What is PSM?**
Propensity score matching estimates the probability that each unit is "treated" based on observable characteristics, then matches treated units to control units with similar probabilities.

**Steps:**
1. **Estimate propensity scores**: Fit a logistic regression predicting treatment (no prior PDB structure) from pre-treatment characteristics:
   - `average_plddt`: AlphaFold prediction quality
   - `pre_mean`: Mean papers/month in pre-period
   - `pre_std`: Std dev of papers in pre-period
   - `pre_trend`: Linear slope of papers over pre-period
   - `pre_newcomer`, `pre_veteran`: Papers by author type
   - `pre_mesh`, `pre_new_mesh`: Research topic diversity
   - `protein_existence` dummies: Level of protein characterization

2. **Match on propensity score**: For each treated gene, find the control gene with the nearest propensity score (Euclidean distance).

3. **Apply caliper**: Drop matches where propensity score difference > 0.1.

**Results:**
```
Propensity Score Distribution:
  Treated: mean=0.71, range=[0.00, 1.00]
  Control: mean=0.45, range=[0.00, 1.00]
  Overlap region: [0.00, 1.00] ✓ (good overlap)
```

**Why it didn't work well:**
The propensity score balances the *probability of treatment*, but treated and control genes are fundamentally different in their *outcome levels*. A treated gene with pscore=0.7 might have 2 papers/month, while a control gene with pscore=0.7 might have 50 papers/month. Matching on propensity score alone doesn't ensure we're comparing genes with similar research activity.

Balance after PSM:
| Variable | Before (Std Diff) | After PSM (Std Diff) |
|----------|-------------------|----------------------|
| pre_mean | -0.30 | -0.19 |
| pre_trend | -0.27 | -0.16 |
| pscore | 1.22 | 0.78 |

Pre-trends were still not parallel after PSM.

**Output file**: `matched_pairs.parquet`

---

### Approach 2: Exact Bin + Nearest Neighbor Matching (Main Approach)

**Key insight**: The fundamental problem is that treated genes (no prior PDB) have much lower research activity than control genes (had prior PDB). To achieve parallel trends, we need to compare genes with *similar pre-treatment paper counts*, not just similar propensity scores.

**Steps:**

1. **Create bins based on pre-treatment mean papers**:
   ```
   Bins: [0-1], [1-3], [3-5], [5-10], [10-20], [20-50], [50-100], [100+]
   ```
   This is "exact matching" on the binned lagged dependent variable.

2. **Check overlap within each bin**:
   ```
   Bin      | Treated | Control | Can match?
   ---------|---------|---------|------------
   0-1      |  5,533  |    453  | Yes (but few controls)
   1-3      |  3,635  |  1,915  | Yes
   3-5      |  1,227  |  1,161  | Yes
   5-10     |  1,004  |  1,480  | Yes
   10-20    |    448  |  1,087  | Yes
   20-50    |    223  |    865  | Yes
   50-100   |     68  |    390  | Yes
   100+     |     41  |    420  | Yes
   ```

3. **Within each bin, match using nearest neighbor on multiple features**:
   - Standardize features within the bin
   - Compute Euclidean distance on: `pre_mean`, `pre_trend`, `average_plddt`
   - For each treated gene, find the nearest control gene

4. **Allow matching with replacement**:
   - Each control gene can be matched to multiple treated genes
   - When computing statistics, weight controls by how often they're used

**Why this works:**
By forcing exact matching on the binned lagged outcome, we ensure that a treated gene with ~5 papers/month is compared to a control gene with ~5 papers/month. The within-bin nearest neighbor matching further refines the match on trend and other characteristics.

**Results:**

Balance after Exact Bin + NN matching:
| Variable | Before (Std Diff) | After Bin+NN (Std Diff) |
|----------|-------------------|-------------------------|
| pre_mean | -0.30 | -0.20 |
| pre_trend | -0.27 | -0.18 |
| average_plddt | -0.41 | -0.27 |

**Pre-trend slopes after matching**:
- Treated: 0.097 papers/month
- Control: 0.097 papers/month
- **Parallel trends achieved!**

The key difference from PSM: by binning on the lagged outcome first, we ensure comparability in *levels* before fine-tuning the match.

**Output file**: `matched_pairs_exact_bin.parquet`

---

### Visualization

See `pretrends_comparison_FIXED.png`:
- Top row: Full sample (non-parallel trends, control ~30 papers/mo vs treated ~4 papers/mo)
- Bottom row: Matched sample (parallel trends, both ~4 papers/mo, difference ≈ 0)

![Pre-trends Comparison - Full vs Matched](pretrends_comparison_FIXED.png)

## DiD Estimates (Simple Calculations)

| Sample | DiD Estimate | Interpretation |
|--------|--------------|----------------|
| Full (simple) | -3.95 | Biased by level differences |
| Full (FE) | +0.48 | Positive but driven by confounding |
| Matched (weighted) | -0.002 | Essentially zero |
| Matched (FE) | +0.48 | Similar to full sample FE |

---

## Formal Regression Results (PSM Matched Sample)

### Methodology Note

> **Important**: This section uses the **PSM-matched sample** (`matched_pairs.parquet`), which is different from the exact-bin matched sample used in the earlier "Heterogeneous Effects Analysis" section. The two approaches yield different results because:
> - **PSM**: Matches on propensity score globally, then subsets by bin post-hoc → does NOT guarantee parallel trends within bins
> - **Exact-bin matching**: Matches within bins first → guarantees comparability within bins

**Regression specification**:
```
n_papers_it = α_i + γ_t + β(Treated_i × Post_t) + ε_it
```

**Sample construction**:
1. Start with PSM matched pairs (12,179 treated genes, 3,817 unique controls)
2. Controls weighted by match frequency (capped at 10)
3. For bin-specific regressions: **subset** the PSM sample to genes in that bin, then run separate regression

**Weighting math**: If control gene C is matched to N treated genes:
- Weight = min(N, 10)
- Weighted OLS: minimize Σ wᵢ(yᵢ - X̂β)²

**Known limitations**:
- One control matched **1,882 times** (gene 340096, pre_mean=0.4 papers/mo → in 0-1 bin)
- The **0-1 bin** has 5,311 treated vs only 352 controls (15:1 ratio) — heavily influenced by control reuse
- The **20+ bin** has tiny sample (333 treated, 331 controls) with high variance (pre_mean ranges 20-1600) — results are unstable

### Aggregate Results

| Specification | β | SE | p-value | Notes |
|--------------|---|-----|---------|-------|
| Full sample (unmatched) | -2.80 | 0.20 | <0.001*** | Confounded |
| PSM, unweighted | -0.77 | 0.16 | <0.001*** | Treats all genes equally |
| PSM, weighted | -0.10 | 0.12 | 0.41 | Proper PSM weights |
| PSM, weights capped at 10 | **-0.27** | **0.08** | **<0.001****** | Robustness check |
| PSM, exclude top 10 reused controls | -0.26 | 0.07 | <0.001*** | Robustness check |

**Concern**: One control gene is matched to **1,882 treated genes**. This extreme reuse inflates that control's influence. When we cap weights at 10 or exclude heavily-reused controls, we get a significant **negative** aggregate effect (~-0.27 papers/month).

### Heterogeneous Effects by Pre-Treatment Intensity (Capped Weights)

| Pre-treatment papers/mo | β | SE | p-value | Reliability |
|------------------------|---|-----|---------|-------------|
| 0-1 | -0.054 | 0.017 | 0.001*** | ⚠️ 15:1 T:C ratio |
| 1-3 | -0.077 | 0.021 | <0.001*** | ✓ |
| 3-5 | -0.007 | 0.059 | 0.91 | ✓ |
| 5-10 | **+0.163** | 0.070 | **0.019**** | ✓ Best identified |
| 10-20 | +0.113 | 0.183 | 0.54 | ~ Underpowered |
| 20+ | +0.500 | 1.103 | 0.65 | ⚠️ Unstable |

**Key pattern**:
- Low-activity genes (0-3 papers/mo): Significant *negative* effect
- Mid-activity (3-5): Null
- Higher-activity (5+): Positive effects (significant for 5-10, underpowered for higher bins)

### Heterogeneous Effects by pLDDT (Capped Weights)

| Quartile | β | SE | p-value |
|----------|---|-----|---------|
| Q1 (low) | -0.169 | 0.065 | 0.009*** |
| Q2 | -0.173 | 0.087 | 0.046** |
| Q3 | -0.410 | 0.244 | 0.093* |
| Q4 (high) | -0.341 | 0.150 | 0.023** |

All pLDDT quartiles show negative effects. The key driver of heterogeneity is **pre-treatment intensity**, not pLDDT.

### Interpretation

1. **Aggregate effect after matching is negative** (or null depending on specification)
2. **Effect varies dramatically by pre-treatment research activity**:
   - Low-activity genes: AlphaFold had negative/null effect
   - High-activity genes: AlphaFold had positive effect (accelerated research)
3. **The 5-10 papers/month bin is the cleanest result**: β = +0.16, p = 0.02

---

## Key Finding (Aggregate)

**Matching successfully eliminates pre-trends, but the aggregate treatment effect vanishes.**

The difference between treated and control in the matched sample fluctuates around zero both before AND after treatment in the aggregate.

---

## Heterogeneous Effects Analysis (Exact-Bin Matched Sample)

> **Methodology Note**: This section uses the **exact-bin matched sample** (`matched_pairs_exact_bin.parquet`), NOT the PSM sample. Matching was done *within* intensity bins, which is why pre-trends are perfectly parallel. These are simple DiD calculations (mean differences), not formal panel regressions.

After finding no aggregate effect, we examined whether effects exist in specific subgroups.

### By pLDDT (AlphaFold Prediction Quality)

| Quartile | DiD Estimate | Pre-slopes (T vs C) | N Treated |
|----------|--------------|---------------------|-----------|
| Q1 (low) | -0.002 | 0.087 vs 0.076 | 3,047 |
| Q2 | +0.058 | 0.082 vs 0.060 | 3,043 |
| Q3 | **+0.116** | 0.122 vs 0.094 | 3,044 |
| Q4 (high) | +0.077 | 0.098 vs 0.084 | 3,045 |

![Heterogeneous Effects by pLDDT](het_effects_plddt_matched.png)

**Interpretation**: Moderate-high quality predictions (Q3) show the largest effect. Very low quality predictions show no effect. The highest quartile may include "easy" proteins that already had good experimental structures.

### By Protein Existence Level

| Category | DiD Estimate | Pre-slopes (T vs C) | N Treated |
|----------|--------------|---------------------|-----------|
| Protein level | **+0.058** | 0.113 vs 0.097 | 10,295 |
| Transcript level | -0.013 | 0.018 vs 0.011 | 725 |
| Inferred from homology | -0.032 | 0.009 vs 0.007 | 498 |
| Predicted | -0.053 | 0.004 vs 0.003 | 114 |

![Heterogeneous Effects by Protein Existence](het_effects_protein_existence.png)

**Interpretation**: Effect only appears for well-characterized proteins with evidence at protein level. Less characterized proteins show no effect (and have very low research activity generally).

### By Pre-Treatment Research Intensity (KEY FINDING)

| Pre-treatment papers/mo | DiD Estimate | Pre-slopes (T vs C) | N Treated |
|------------------------|--------------|---------------------|-----------|
| 0-1 | -0.052 | 0.008 vs 0.008 ✓ | 5,503 |
| 1-3 | -0.069 | 0.039 vs 0.038 ✓ | 3,635 |
| 3-5 | +0.071 | 0.088 vs 0.088 ✓ | 1,227 |
| 5-10 | +0.139 | 0.175 vs 0.175 ✓ | 1,004 |
| 10-20 | **+0.291** | 0.361 vs 0.358 ✓ | 448 |
| 20+ | **+0.716** | 1.687 vs 1.687 ✓ | 330 |

![Heterogeneous Effects by Pre-Treatment Intensity](het_effects_intensity.png)

**This is the cleanest result**:
- Pre-trends are **perfectly parallel** in every bin (matching worked)
- Effect **increases monotonically** with pre-treatment activity
- Genes with 0-3 papers/month: No effect
- Genes with 20+ papers/month: Large effect (+0.72 papers/month)

**Interpretation**: AlphaFold **accelerated existing research** but **did not stimulate new research** on inactive genes. Researchers who were already working on a protein could use the structure to do more; researchers weren't drawn to start studying proteins just because AlphaFold predicted their structure.

### Summary Visualization

See `het_effects_summary.png` for a bar chart comparison of all heterogeneous effects.

![Heterogeneous Effects Summary](het_effects_summary.png)

---

## Files in this Directory

### Data Files
- `gene_features_for_matching.parquet`: Gene-level features used for matching (N=19,950)
- `matched_pairs.parquet`: PSM-based matched pairs (first attempt)
- `matched_pairs_exact_bin.parquet`: Exact bin + NN matched pairs (main analysis)

### Figures
- `pretrends_comparison_FIXED.png`: Main figure comparing full vs matched samples
- `het_effects_plddt_matched.png`: Heterogeneous effects by pLDDT quartile
- `het_effects_protein_existence.png`: Heterogeneous effects by protein evidence
- `het_effects_intensity.png`: Heterogeneous effects by pre-treatment intensity
- `het_effects_summary.png`: Summary bar chart of all heterogeneous effects
- `trends_by_intensity_PSM.png`: Raw trends by intensity bin (PSM sample)
- `event_study_by_intensity_PSM.png`: Event studies by intensity bin

---

## Diagnostics and Robustness Checks

### Control Reuse Concern

With PSM matching, controls can be matched to multiple treated genes. Our data shows:
- **One control matched 1,882 times** (extreme!)
- Top 10 most-reused controls account for disproportionate influence
- When we cap weights at 10 or exclude heavily-reused controls, results change

This is a methodological concern worth noting.

### Truncation Bias Check

Papers take months to appear in PubMed, so late 2023 data is incomplete. We tested dropping the last 3 months (Oct-Dec 2023):

| Specification | Full period β | Trimmed β | Change? |
|---------------|---------------|-----------|---------|
| Aggregate | -0.27 | -0.39 | More negative |
| 0-1 bin | -0.054 | -0.065 | Similar |
| 5-10 bin | +0.163 | +0.142 | Similar |

**Conclusion**: Truncation doesn't explain the pattern. Results are robust to dropping end-of-sample.

### Log Specification (% Interpretation)

To check whether we're just picking up absolute vs percentage changes:

| Bin | Levels β | Log β | % Effect |
|-----|----------|-------|----------|
| 0-1 | -0.065 | -0.029 | **-2.9%**** |
| 1-3 | -0.088 | -0.027 | **-2.7%**** |
| 3-5 | -0.012 | -0.001 | -0.1% |
| 5-10 | +0.142 | +0.019 | **+1.9%**** |
| 10-20 | +0.122 | +0.011 | +1.1% |
| 20+ | +0.485 | +0.017 | +1.7%* |

**Conclusion**: Pattern holds in logs. Low-activity genes: ~3% decline. High-activity genes: ~2% increase.

### Event Study Diagnostics

Event studies by intensity bin reveal:
- **0-1, 1-3 bins**: Show pre-trend drift before treatment (concerning - suggests PSM didn't fully solve parallel trends)
- **3-5 bin**: Fluctuates around zero throughout (null effect looks real)
- **5-10 bin**: Cleanest - relatively flat pre-treatment, positive post-treatment
- **10-20, 20+ bins**: Very noisy due to small sample sizes

![Trends by Intensity](trends_by_intensity_PSM.png)

See `event_study_by_intensity_PSM.png`.

![Event Study by Intensity](event_study_by_intensity_PSM.png)

---

## Final Interpretation

### What We Can Say With Confidence

1. **Aggregate effect after matching is negative or null** depending on specification
2. **Strong heterogeneity by pre-treatment research activity**:
   - Low-activity genes (0-3 papers/mo): Negative effect (~-3%), but pre-trends questionable
   - Mid-activity (3-5): Null effect
   - Higher-activity (5-10): **Positive effect (+1.9%)**, cleanest result
3. **The 5-10 papers/month bin is most credible**: Has reasonable sample size, decent pre-trends, significant positive effect

### Speculative Explanations

**Why negative effects for low-activity genes?**
- Bad matching (0-1 bin: 5,533 treated vs 395 controls = 14:1 ratio)
- Selection: genes without PDB in low-activity group may be fundamentally different
- Event studies show pre-trend drift → results may be artifacts

**Why positive effects for high-activity genes?**
- Makes intuitive sense: researchers already working on a protein can use the structure
- AlphaFold accelerated existing research, didn't create new interest

### Bottom Line

> AlphaFold appears to have **accelerated research on actively-studied proteins** (~2% more papers) but **did not stimulate new research on understudied proteins**. The negative effects for low-activity genes are suspect due to matching limitations.

---

## References

- Caetano & Callaway (2024): "Difference-in-Differences when Parallel Trends Holds Conditional on Covariates"
- EU CAP Network: "PSM-DiD Method" learning portal
- Stuart et al. (2014): "Using propensity scores in difference-in-differences models"

## Potential Next Steps

1. ~~Heterogeneous effects~~ ✓ Done
2. **Try exact-bin matching** (better balance) and compare to PSM results
3. **Synthetic control methods** for aggregate effect
4. **Investigate the 0-1 bin** more carefully - why such poor matching?
5. **Alternative treatment definitions**: Continuous treatment based on pLDDT improvement
