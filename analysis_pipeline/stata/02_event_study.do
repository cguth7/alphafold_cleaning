/*******************************************************************************
02_event_study.do - Event Study Regressions
================================================================================

Runs weighted event study regressions using the CEM-matched panel.

SPECIFICATIONS:
  1. LEVELS:  Y ~ treated×period FE + gene FE + semester FE
  2. DY:      ΔY ~ treated×period FE + gene FE + semester FE  [PREFERRED]
  3. DASINH:  Δasinh(Y) ~ treated×period FE + gene FE + semester FE

OUTPUT:
  - CSV files with coefficients for each outcome × specification
  - Event study figures (PNG)

*******************************************************************************/

di as txt ""
di as txt "============================================================"
di as txt "STEP 2: EVENT STUDY REGRESSIONS"
di as txt "============================================================"

*==============================================================================*
* Load matched panel
*==============================================================================*

use "${OUTPUT_PATH}/derived/matched_panel_semester.dta", clear

* Get the factor-variable base period
qui sum rel_semester, meanonly
local min_sem = r(min)
local base_sem_fv = ${BASE_SEM} - `min_sem'

di as txt "  Base period (original): ${BASE_SEM}"
di as txt "  Base period (FV-safe):  `base_sem_fv'"
di as txt ""

*==============================================================================*
* Event study program
*==============================================================================*

capture program drop run_eventstudy
program define run_eventstudy
    syntax varname, spec(string) outcome(string) base(integer) shift(integer)

    * Run regression with two-way FE
    qui reghdfe `varlist' ib`base'.rel_sem_fv##i.treated [aw=cem_weight], ///
        absorb(gene_id semester) vce(cluster gene_id)

    * Get coefficient count
    qui levelsof rel_sem_fv if e(sample), local(periods)

    * Create results file
    tempfile results
    postfile handle int period double coef se ci_low ci_high using `results', replace

    foreach k of local periods {
        if `k' == `base' {
            post handle (`k') (0) (0) (0) (0)
        }
        else {
            local b = _b[1.treated#`k'.rel_sem_fv]
            local se = _se[1.treated#`k'.rel_sem_fv]
            local lo = `b' - 1.96 * `se'
            local hi = `b' + 1.96 * `se'
            post handle (`k') (`b') (`se') (`lo') (`hi')
        }
    }
    postclose handle

    * Load and adjust period
    preserve
        use `results', clear
        replace period = period + `shift'
        sort period

        * Save CSV
        export delimited using "${OUTPUT_PATH}/eventstudy_`spec'_`outcome'.csv", replace

        * Create figure
        twoway (rcap ci_high ci_low period, lcolor(navy)) ///
               (connected coef period, mcolor(navy) lcolor(navy)), ///
               yline(0, lcolor(black) lpattern(solid)) ///
               xline(${BASE_SEM}, lcolor(red) lpattern(dash)) ///
               title("`outcome' - `spec' (Semester)") ///
               xtitle("Semesters relative to AlphaFold2") ///
               ytitle("Treatment Effect") ///
               legend(off) ///
               graphregion(color(white)) bgcolor(white)

        graph export "${FIGURE_PATH}/eventstudy_`spec'_`outcome'.png", replace width(1200)
    restore
end

*==============================================================================*
* Run event studies for each outcome
*==============================================================================*

qui sum rel_semester, meanonly
local shift = r(min)

foreach y of global OUTCOMES {

    di as txt "  Outcome: `y'"

    * Check if variable exists
    cap confirm variable `y'
    if _rc {
        di as txt "    Skipping `y' (not found)"
        continue
    }

    *--- LEVELS specification ---*
    di as txt "    Running LEVELS..."
    run_eventstudy `y', spec("LEVELS") outcome("`y'") base(`base_sem_fv') shift(`shift')

    *--- DY specification ---*
    cap confirm variable D_`y'
    if !_rc {
        di as txt "    Running DY..."
        run_eventstudy D_`y', spec("DY") outcome("`y'") base(`base_sem_fv') shift(`shift')
    }

    *--- DASINH specification ---*
    cap confirm variable D_asinh_`y'
    if !_rc {
        di as txt "    Running DASINH..."
        run_eventstudy D_asinh_`y', spec("DASINH") outcome("`y'") base(`base_sem_fv') shift(`shift')
    }

    di as txt ""
}

di as txt ">>> Event studies complete."
di as txt "    CSV files: ${OUTPUT_PATH}/eventstudy_*.csv"
di as txt "    Figures:   ${FIGURE_PATH}/eventstudy_*.png"
di as txt ""
