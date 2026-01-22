/*******************************************************************************
01_cem_matching.do - Coarsened Exact Matching for AlphaFold2 Analysis
================================================================================

Creates matched treatment/control samples using CEM.

TREATMENT DEFINITION:
  treated = 1 if num_deposits == 0 (genes with NO prior PDB structures)

MATCHING VARIABLES (pre-treatment only):
  - pre_mean_papers:   mean publications before July 2021
  - pre_mean_newcom:   mean newcomer papers
  - pre_mean_veteran:  mean veteran papers
  - pre_mean_top10:    mean top-10% papers
  - pre_sd_papers:     std dev of papers

NOTE: pLDDT is NOT used for matching (per Matteo - it's AF2-generated)

OUTPUT:
  - derived/matched_panel_semester.dta (main analysis file)
  - derived/gene_level_features.dta (gene-level summaries)

*******************************************************************************/

di as txt ""
di as txt "============================================================"
di as txt "STEP 1: CEM MATCHING"
di as txt "============================================================"

*==============================================================================*
* 1. Load data and create time variables
*==============================================================================*

di as txt "Loading data..."
use "${DATA_PATH}", clear

* Get min ym for sequential indexing
qui sum ym, meanonly
local ym_min = r(min)

* Sequential month index
gen long ym_seq = ym - `ym_min' + 1
label var ym_seq "Sequential month index (1..T)"

* Treatment sequence position
local treatment_seq = ${TREATMENT_YM} - `ym_min' + 1

* Treatment indicator
gen byte treated = (num_deposits == 0)
label var treated "Treatment group (no prior PDB deposits)"
label define treated_lbl 0 "Control (has PDB)" 1 "Treated (no PDB)"
label values treated treated_lbl

* Post-treatment indicator
gen byte post = (ym_seq >= `treatment_seq')
label var post "Post-treatment period"

* Quarter and semester indices
gen int quarter  = floor((ym_seq - 1) / 3)
gen int semester = floor((ym_seq - 1) / 6)
label var quarter  "Quarter index from start"
label var semester "Semester index from start"

* Relative time (event time)
local treatment_quarter  = floor((`treatment_seq' - 1) / 3)
local treatment_semester = floor((`treatment_seq' - 1) / 6)

gen int rel_quarter  = quarter  - `treatment_quarter'
gen int rel_semester = semester - `treatment_semester'
label var rel_quarter  "Event time in quarters (0 = treatment)"
label var rel_semester "Event time in semesters (0 = treatment)"

di as txt "  Treatment month (ym_seq): `treatment_seq'"
di as txt "  Treatment semester: `treatment_semester'"

* Save base panel
save "${OUTPUT_PATH}/derived/panel_base.dta", replace

*==============================================================================*
* 2. Compute pre-treatment gene-level features
*==============================================================================*

di as txt "Computing pre-treatment features..."

preserve
    * Keep only pre-treatment periods
    keep if ym_seq < `treatment_seq'

    * Collapse to gene level
    collapse (mean) pre_mean_papers = n_papers ///
             (mean) pre_mean_newcom = n_newcomer_papers ///
             (mean) pre_mean_veteran = n_veteran_papers ///
             (mean) pre_mean_top10 = n_top10_y ///
             (sd)   pre_sd_papers = n_papers, ///
             by(gene_id)

    * Fill missing SD with 0
    replace pre_sd_papers = 0 if missing(pre_sd_papers)

    * Label variables
    label var pre_mean_papers  "Pre-treatment mean: total papers"
    label var pre_mean_newcom  "Pre-treatment mean: newcomer papers"
    label var pre_mean_veteran "Pre-treatment mean: veteran papers"
    label var pre_mean_top10   "Pre-treatment mean: top-10% papers"
    label var pre_sd_papers    "Pre-treatment SD: total papers"

    save "${OUTPUT_PATH}/derived/pre_features.dta", replace
restore

* Merge pre-features back
merge m:1 gene_id using "${OUTPUT_PATH}/derived/pre_features.dta", nogen

* Get gene-level treatment status and other time-invariant vars
preserve
    collapse (first) treated average_plddt gene_name protein_id, by(gene_id)
    merge 1:1 gene_id using "${OUTPUT_PATH}/derived/pre_features.dta", nogen

    * Fill missing
    foreach v in pre_mean_papers pre_mean_newcom pre_mean_veteran pre_mean_top10 pre_sd_papers {
        replace `v' = 0 if missing(`v')
    }
    replace average_plddt = 0 if missing(average_plddt)

    save "${OUTPUT_PATH}/derived/gene_level_features.dta", replace
restore

*==============================================================================*
* 3. Run CEM matching
*==============================================================================*

di as txt "Running CEM matching..."

use "${OUTPUT_PATH}/derived/gene_level_features.dta", clear

* Coarsen variables into quantile bins
foreach v in pre_mean_papers pre_mean_newcom pre_mean_veteran pre_mean_top10 pre_sd_papers {
    xtile cem_`v' = `v', nq(${CEM_NQ})
}

* Run CEM (NOTE: pLDDT intentionally excluded per Matteo's instruction)
cem cem_pre_mean_papers cem_pre_mean_newcom cem_pre_mean_veteran ///
    cem_pre_mean_top10 cem_pre_sd_papers, treatment(treated)

* Check for matched indicator
cap confirm variable cem_matched
if _rc {
    di as error "CEM did not create cem_matched variable. Check cem output."
    exit 459
}

* Keep matched genes only
keep if cem_matched == 1

* Rename weight
gen double cem_weight = cem_weights
label var cem_weight "CEM matching weight"

* Report sample sizes
count if treated == 1
local n_treated = r(N)
count if treated == 0
local n_control = r(N)

di as txt ""
di as txt "  CEM Matching Results:"
di as txt "    Treated genes: `n_treated'"
di as txt "    Control genes: `n_control'"
di as txt "    Total matched: `=`n_treated' + `n_control''"

* Save matched genes
keep gene_id treated cem_weight cem_strata
save "${OUTPUT_PATH}/derived/matched_genes.dta", replace

*==============================================================================*
* 4. Create matched semester panel
*==============================================================================*

di as txt "Creating matched semester panel..."

* Load base panel
use "${OUTPUT_PATH}/derived/panel_base.dta", clear

* Keep only matched genes
merge m:1 gene_id using "${OUTPUT_PATH}/derived/matched_genes.dta", keep(match) nogen

* Collapse to semester level
collapse (sum) n_papers n_newcomer_papers n_veteran_papers n_top10_y ///
         (first) treated cem_weight gene_name protein_id average_plddt ///
                 num_deposits rel_semester post, ///
         by(gene_id semester)

* Create delta (first difference) variables
sort gene_id semester
foreach y of global OUTCOMES {
    by gene_id: gen double D_`y' = `y' - `y'[_n-1]
    label var D_`y' "First difference: `y'"

    gen double asinh_`y' = asinh(`y')
    by gene_id: gen double D_asinh_`y' = asinh_`y' - asinh_`y'[_n-1]
    label var asinh_`y' "Asinh transform: `y'"
    label var D_asinh_`y' "Delta asinh: `y'"
}

* Create factor-variable safe relative semester (non-negative)
qui sum rel_semester, meanonly
local min_sem = r(min)
gen int rel_sem_fv = rel_semester - `min_sem'
label var rel_sem_fv "Relative semester (factor-variable safe)"

* Linear time trend
gen double t_sem = rel_semester
label var t_sem "Linear time trend"

* Pre-period indicator
gen byte pre = (rel_semester < 0)
label var pre "Pre-treatment period"

* Report panel size
qui count
local n_obs = r(N)
qui unique gene_id
local n_genes = r(unique)

di as txt ""
di as txt "  Semester Panel Created:"
di as txt "    Observations: `n_obs'"
di as txt "    Genes: `n_genes'"
di as txt "    Semesters: `=`n_obs'/`n_genes''"

* Save
save "${OUTPUT_PATH}/derived/matched_panel_semester.dta", replace

di as txt ""
di as txt "  Saved: ${OUTPUT_PATH}/derived/matched_panel_semester.dta"
di as txt ""
di as txt ">>> CEM Matching complete."
di as txt ""
