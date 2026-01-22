# Honest Assessment: Is This Analysis Meaningful?

## TL;DR

**Your initial skepticism was justified.** The analysis has serious problems that make the results hard to interpret causally. The "flat pre-trends" in the ΔY specification are somewhat illusory - they're flat because both groups were growing, and the "positive effect" post-treatment is partly driven by both groups declining (with control declining more).

---

## The Core Problems

### 1. Massive Baseline Imbalance (CEM Didn't Fix It)

| Group | Papers/Semester | Description |
|-------|-----------------|-------------|
| Control | ~180 | Well-studied genes (have PDB structures) |
| Treated | ~25 | Obscure genes (no PDB structures) |

**The 7x difference persists after CEM matching.** The matching balanced on pre-treatment *trends* (means, SDs), but control genes are fundamentally different - they're the "popular" genes that big labs study.

### 2. The "Effect" Is Actually About Decline, Not Growth

Raw data shows:
```
                    Sem -3   Sem -1   Sem 2    Sem 4
Treated genes:        20       28       29       18
Control genes:       143      204      217      142
```

**Both groups peaked around Semester 2 and then DECLINED sharply.**

The ΔY specification shows a "positive effect" because:
- Both groups are declining
- Control genes decline MORE (from 217 → 142 = -35%)
- Treated genes decline LESS (from 29 → 18 = -38%)
- Wait, actually treated declined slightly more in percentage terms!

The "effect" comes from the ABSOLUTE decline being larger for control (75 papers vs 11 papers), which shows up as a positive coefficient when we difference.

### 3. Data Truncation in 2023

Semester 4 (Jul-Dec 2023) shows a ~35% drop for BOTH groups. This is almost certainly **incomplete data** - papers published in late 2023 haven't made it into PubTator yet.

This truncation is DRIVING the "positive effects" in periods 3-4.

### 4. COVID Confounds Everything

The pre-treatment period (Semesters -3 to -1) spans:
- Jan 2020: Pre-COVID
- Jul 2020: COVID peak
- Jan 2021: Vaccination begins

Both groups were growing during this period, but for reasons unrelated to AF2.

### 5. Only 8 Semesters Total

- 3 pre-treatment periods
- 5 post-treatment periods
- That's barely enough to estimate anything with confidence

---

## Why Danilo's Approach "Works"

Danilo's key insight was that ΔY + CEM matching produces "flat pre-trends." But this is partly mechanical:

1. **ΔY removes level differences** - if treated has 25 papers and control has 180, differencing removes that gap
2. **Both groups were trending up pre-treatment** - so the difference-in-differences is ~0
3. **Both groups trend down post-treatment** - but control drops more in absolute terms → positive DiD

This isn't necessarily wrong, but the interpretation "AF2 caused more publications" is misleading. A more accurate interpretation is "AF2 genes experienced less decline than control genes in 2022-2023."

---

## Stata vs Python Comparison

| Aspect | Danilo (Stata) | Our (Python) | Verdict |
|--------|----------------|--------------|---------|
| Method | reghdfe (proper 2-way FE) | Simple DiD | Different |
| Coefficient magnitudes | Small (0.01-0.1) | Large (1-3) | 20-100x different |
| Qualitative pattern | Flat pre, positive post | Flat pre, positive post | Same |
| Statistical significance | Tight SEs | Wider SEs | His is tighter |

**The magnitudes differ because reghdfe properly absorbs gene and time fixed effects.** Our simple DiD doesn't account for the massive variation in gene-level publication rates.

---

## What the DOI Data Can Do

**Ready for Danilo:**
- `analysis_pipeline/data/dois_for_shi_evans.csv`
- 1.37 million unique DOIs (2020-2023)
- He can run his Shi & Evans model on these

**We can't do more without:**
- His trained novelty model (via collaborator Filippo)
- Or building our own novelty measure (would need OpenAlex API for abstracts/citations)

---

## What Would Make This Analysis More Convincing?

1. **Better control group**: Match on gene "importance" or research intensity, not just pre-trends
2. **Longer pre-period**: Need 5+ pre-treatment periods to really test parallel trends
3. **Account for data truncation**: Either update data or drop 2023
4. **COVID controls**: Add COVID-era indicators or restrict to post-vaccine period
5. **Different treatment definition**: Maybe use AF2 prediction quality (pLDDT) rather than PDB deposits
6. **Placebo tests**: Run the same analysis with fake treatment dates

---

## Bottom Line for Your Meeting

**What to tell Matteo:**

1. "I replicated Danilo's analysis in Python - the qualitative patterns match"

2. "But when I dug into the numbers, I have concerns:
   - The 'effect' is driven by both groups declining, not treated growing
   - There's likely data truncation in 2023
   - CEM didn't fix the 7x baseline imbalance
   - We only have 8 semesters of data"

3. "The DOI file is ready for Danilo's Shi & Evans analysis"

4. "I think we need to discuss whether this is publishable or if we're p-hacking our way to flat pre-trends"

**Your original instinct was right** - this feels like "trying stuff until it works" rather than a pre-specified analysis. The ΔY specification is defensible in principle, but the combination of matching + differencing + short time window + data truncation makes the results hard to interpret causally.
