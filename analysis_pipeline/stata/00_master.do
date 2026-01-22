/*******************************************************************************
ALPHAFOLD2 IMPACT ANALYSIS - MASTER PIPELINE
================================================================================

This master file runs the complete Stata analysis pipeline:
  1. CEM matching and panel construction
  2. Event study regressions
  3. Diagnostics and honest assessment

Based on Danilo Messinese's original code, cleaned up with added diagnostics.

USAGE:
  1. Set the PROJECT_PATH global below to your project root
  2. Run this file: do "analysis_pipeline/stata/00_master.do"

REQUIREMENTS:
  - Stata 15+ (for reghdfe)
  - Packages: cem, reghdfe, ftools, estout (auto-installed if missing)

*******************************************************************************/

clear all
set more off
set varabbrev off

*==============================================================================*
* CONFIGURATION - EDIT THIS SECTION
*==============================================================================*

* Project root (where final_panel_CLEAN.dta lives in final/ subfolder)
* CHANGE THIS TO YOUR PATH:
global PROJECT_PATH "/Users/maxguthmann/Downloads/Development/Work/Alphafold_2"

* Treatment timing
global TREATMENT_YM   738    // July 2021 in Stata ym format

* Matching parameters
global N_STRATA       20     // Bins for stratification
global CEM_NQ         10     // Quantile bins for CEM

* Outcomes to analyze
global OUTCOMES "n_papers n_newcomer_papers n_veteran_papers n_top10_y"

* Base period for event study (relative semester)
global BASE_SEM -1

*==============================================================================*
* DERIVED PATHS (don't edit)
*==============================================================================*

global DATA_PATH      "${PROJECT_PATH}/final/final_panel_CLEAN.dta"
global STATA_PATH     "${PROJECT_PATH}/analysis_pipeline/stata"
global OUTPUT_PATH    "${PROJECT_PATH}/analysis_pipeline/stata/output"
global FIGURE_PATH    "${PROJECT_PATH}/analysis_pipeline/stata/figures"

* Create output directories
cap mkdir "${OUTPUT_PATH}"
cap mkdir "${FIGURE_PATH}"
cap mkdir "${OUTPUT_PATH}/derived"

*==============================================================================*
* INSTALL REQUIRED PACKAGES
*==============================================================================*

foreach pkg in cem reghdfe ftools estout unique {
    cap which `pkg'
    if _rc {
        di as txt "Installing `pkg'..."
        ssc install `pkg', replace
    }
}

* reghdfe needs ftools
cap reghdfe, compile

*==============================================================================*
* RUN PIPELINE
*==============================================================================*

di as txt ""
di as txt "============================================================"
di as txt "ALPHAFOLD2 IMPACT ANALYSIS PIPELINE"
di as txt "============================================================"
di as txt "Project path: ${PROJECT_PATH}"
di as txt "Data file:    ${DATA_PATH}"
di as txt "============================================================"
di as txt ""

* Step 1: CEM Matching
di as txt ">>> STEP 1: CEM Matching"
do "${STATA_PATH}/01_cem_matching.do"

* Step 2: Event Studies
di as txt ">>> STEP 2: Event Study Regressions"
do "${STATA_PATH}/02_event_study.do"

* Step 3: Diagnostics
di as txt ">>> STEP 3: Diagnostics & Honest Assessment"
do "${STATA_PATH}/03_diagnostics.do"

di as txt ""
di as txt "============================================================"
di as txt "PIPELINE COMPLETE"
di as txt "============================================================"
di as txt "Outputs saved to: ${OUTPUT_PATH}"
di as txt "Figures saved to: ${FIGURE_PATH}"
di as txt "============================================================"
