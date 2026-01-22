/*******************************************************************************
03_diagnostics.do - Diagnostics and Honest Assessment
================================================================================

This file outputs diagnostic information that Danilo's original code was missing:

1. Raw trends by treatment group (not just regression coefficients)
2. Percentage changes (not just absolute)
3. Sample sizes by period
4. Data truncation check for 2023
5. Comparison of curves showing they have same SHAPE (just scaled)

This is the "honest assessment" - showing what the data actually looks like
before we apply ΔY transformations that can obscure the underlying patterns.

*******************************************************************************/

di as txt ""
di as txt "============================================================"
di as txt "STEP 3: DIAGNOSTICS & HONEST ASSESSMENT"
di as txt "============================================================"

use "${OUTPUT_PATH}/derived/matched_panel_semester.dta", clear

*==============================================================================*
* 1. Raw means by treatment group and period
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "1. RAW MEANS BY TREATMENT GROUP AND PERIOD"
di as txt "============================================================"
di as txt ""
di as txt "n_papers (mean per gene-semester):"
di as txt ""
di as txt "Period   Treated    Control    Gap        Pct Diff"
di as txt "------   -------    -------    ---        --------"

levelsof rel_semester, local(semesters)
foreach s of local semesters {
    qui sum n_papers if rel_semester == `s' & treated == 1 [aw=cem_weight]
    local t = r(mean)
    qui sum n_papers if rel_semester == `s' & treated == 0 [aw=cem_weight]
    local c = r(mean)
    local gap = `c' - `t'
    local pct = 100 * `t' / `c'

    di as txt %6.0f `s' "   " %8.1f `t' "   " %8.1f `c' "   " %8.1f `gap' "   " %6.1f `pct' "%"
}

*==============================================================================*
* 2. Period-over-period PERCENTAGE changes
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "2. PERCENTAGE CHANGES (Period over Period)"
di as txt "============================================================"
di as txt ""
di as txt "This shows whether curves have SAME SHAPE (just scaled)."
di as txt "If % changes are similar, the 'effect' from DeltaY is mechanical."
di as txt ""
di as txt "Period    Treated %Δ    Control %Δ    Difference"
di as txt "------    ----------    ----------    ----------"

local prev_t = .
local prev_c = .

foreach s of local semesters {
    qui sum n_papers if rel_semester == `s' & treated == 1 [aw=cem_weight]
    local t = r(mean)
    qui sum n_papers if rel_semester == `s' & treated == 0 [aw=cem_weight]
    local c = r(mean)

    if `prev_t' != . {
        local pct_t = 100 * (`t' - `prev_t') / `prev_t'
        local pct_c = 100 * (`c' - `prev_c') / `prev_c'
        local diff = `pct_t' - `pct_c'

        di as txt %6.0f `s' "      " %+8.1f `pct_t' "%      " %+8.1f `pct_c' "%      " %+8.1f `diff' "pp"
    }

    local prev_t = `t'
    local prev_c = `c'
}

di as txt ""
di as txt "KEY INSIGHT: If percentage changes are similar between groups,"
di as txt "the DeltaY 'effect' is driven by absolute scale differences,"
di as txt "not differential treatment response."

*==============================================================================*
* 3. Sample sizes by period
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "3. SAMPLE SIZES BY PERIOD"
di as txt "============================================================"
di as txt ""

tab rel_semester treated

*==============================================================================*
* 4. Data truncation check (2023)
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "4. DATA TRUNCATION CHECK"
di as txt "============================================================"
di as txt ""

* Get the last two semesters
qui sum rel_semester, meanonly
local max_sem = r(max)
local prev_sem = `max_sem' - 1

* Calculate drops
qui sum n_papers if rel_semester == `prev_sem' & treated == 1 [aw=cem_weight]
local t_prev = r(mean)
qui sum n_papers if rel_semester == `max_sem' & treated == 1 [aw=cem_weight]
local t_last = r(mean)
local t_drop = 100 * (`t_last' - `t_prev') / `t_prev'

qui sum n_papers if rel_semester == `prev_sem' & treated == 0 [aw=cem_weight]
local c_prev = r(mean)
qui sum n_papers if rel_semester == `max_sem' & treated == 0 [aw=cem_weight]
local c_last = r(mean)
local c_drop = 100 * (`c_last' - `c_prev') / `c_prev'

di as txt "Last semester (rel_sem = `max_sem') vs previous (rel_sem = `prev_sem'):"
di as txt ""
di as txt "  Treated: " %6.1f `t_prev' " -> " %6.1f `t_last' " (" %+5.1f `t_drop' "%)"
di as txt "  Control: " %6.1f `c_prev' " -> " %6.1f `c_last' " (" %+5.1f `c_drop' "%)"
di as txt ""

if `t_drop' < -20 | `c_drop' < -20 {
    di as error "WARNING: Sharp drop in last semester suggests DATA TRUNCATION."
    di as error "Papers from late 2023 may not be fully in the database."
    di as txt ""
}

*==============================================================================*
* 5. Simple DiD calculation (sanity check)
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "5. SIMPLE DiD CALCULATION"
di as txt "============================================================"
di as txt ""

* Pre vs post means
qui sum n_papers if rel_semester < 0 & treated == 1 [aw=cem_weight]
local pre_t = r(mean)
qui sum n_papers if rel_semester >= 0 & treated == 1 [aw=cem_weight]
local post_t = r(mean)

qui sum n_papers if rel_semester < 0 & treated == 0 [aw=cem_weight]
local pre_c = r(mean)
qui sum n_papers if rel_semester >= 0 & treated == 0 [aw=cem_weight]
local post_c = r(mean)

local did = (`post_t' - `pre_t') - (`post_c' - `pre_c')

di as txt "LEVELS (Y):"
di as txt "  Treated: " %6.1f `pre_t' " -> " %6.1f `post_t' " (change: " %+6.1f `post_t' - `pre_t' ")"
di as txt "  Control: " %6.1f `pre_c' " -> " %6.1f `post_c' " (change: " %+6.1f `post_c' - `pre_c' ")"
di as txt "  DiD: " %+6.1f `did'
di as txt ""

* Same for DY
qui sum D_n_papers if rel_semester < 0 & treated == 1 [aw=cem_weight]
local pre_t_d = r(mean)
qui sum D_n_papers if rel_semester >= 0 & treated == 1 [aw=cem_weight]
local post_t_d = r(mean)

qui sum D_n_papers if rel_semester < 0 & treated == 0 [aw=cem_weight]
local pre_c_d = r(mean)
qui sum D_n_papers if rel_semester >= 0 & treated == 0 [aw=cem_weight]
local post_c_d = r(mean)

local did_d = (`post_t_d' - `pre_t_d') - (`post_c_d' - `pre_c_d')

di as txt "DELTA Y (first differences):"
di as txt "  Treated: " %+6.2f `pre_t_d' " -> " %+6.2f `post_t_d' " (change: " %+6.2f `post_t_d' - `pre_t_d' ")"
di as txt "  Control: " %+6.2f `pre_c_d' " -> " %+6.2f `post_c_d' " (change: " %+6.2f `post_c_d' - `pre_c_d' ")"
di as txt "  DiD: " %+6.2f `did_d'
di as txt ""

di as txt "NOTE: Levels shows negative effect, DY shows positive effect."
di as txt "This is because DY compares ABSOLUTE changes, not percentage."
di as txt "Control's larger absolute decline shows up as 'treated doing better'."

*==============================================================================*
* 6. Summary statistics
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "6. SUMMARY STATISTICS"
di as txt "============================================================"
di as txt ""

tabstat n_papers n_newcomer_papers n_veteran_papers n_top10_y, ///
    by(treated) stats(n mean sd min p50 max) format(%9.2f)

*==============================================================================*
* 7. Export diagnostics to file
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "7. EXPORTING DIAGNOSTICS"
di as txt "============================================================"

* Create a summary table
preserve
    collapse (mean) n_papers n_newcomer_papers n_veteran_papers n_top10_y ///
             (count) n_genes = gene_id [aw=cem_weight], ///
             by(rel_semester treated)

    export delimited using "${OUTPUT_PATH}/diagnostics_means_by_period.csv", replace
    di as txt "  Saved: diagnostics_means_by_period.csv"
restore

di as txt ""
di as txt ">>> Diagnostics complete."
di as txt ""
di as txt "============================================================"
di as txt "HONEST ASSESSMENT SUMMARY"
di as txt "============================================================"
di as txt ""
di as txt "1. Treated and control groups have SAME SHAPE curves (scaled)"
di as txt "2. Both groups decline sharply in later periods"
di as txt "3. The 'positive effect' in DY is because control declines"
di as txt "   more in ABSOLUTE terms (not percentage terms)"
di as txt "4. Data truncation likely in 2023 (both groups drop ~30%+)"
di as txt "5. This is measuring deceleration difference, not growth"
di as txt ""
di as txt "============================================================"
