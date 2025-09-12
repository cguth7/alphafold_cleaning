* ==============================================================================
* AlphaFold Impact Analysis Template
* Updated for final_panel_COMPLETE dataset structure
* ==============================================================================

* Install required packages
ssc install estout
ssc install ftools
ssc install reghdfe
ssc install regsave

* Load the final balanced panel
* Note: Generate .dta file first if needed:
* python: df = pd.read_parquet('final/final_panel_COMPLETE.parquet'); df.to_stata('final/final_panel_COMPLETE.dta', write_index=False)
use "final/final_panel_COMPLETE.dta", clear

* Check dataset structure
describe
summarize ym year month

* ==============================================================================
* DATA PREPARATION
* ==============================================================================

* Create quarterly time variables for analysis
* ym is in YYYYMM format (202001-202312), need to convert to Stata quarterly
gen year_temp = int(ym/100)
gen month_temp = ym - year_temp*100

* Create proper Stata quarterly date
gen quarter_temp = quarter(mdy(month_temp, 1, year_temp))
gen quarter_num = yq(year_temp, quarter_temp)
format quarter_num %tq

* Clean up temporary variables
drop year_temp month_temp quarter_temp

* Verify quarterly conversion worked
tab quarter_num, missing
list ym year month quarter_num in 1/20

* ==============================================================================
* ALPHAFOLD TREATMENT SETUP
* ==============================================================================

* Define treatment based on AlphaFold confidence
* High confidence proteins are "treated" by AlphaFold
gen high_confidence = (average_plddt > 70) if !missing(average_plddt)
replace high_confidence = 0 if missing(high_confidence)

* Alternative: Use protein existence levels as treatment
* gen high_evidence = (protein_existence == "Evidence at protein level") if !missing(protein_existence)

* Define post-AlphaFold periods
* AlphaFold DB launched July 2021 ≈ 2021q3
* Major expansion July 2022 ≈ 2022q3
gen after_af_launch = (quarter_num >= tq(2021q3))
gen after_af_expansion = (quarter_num >= tq(2022q3))

* Treatment interactions
gen treat_launch = high_confidence * after_af_launch
gen treat_expansion = high_confidence * after_af_expansion

* ==============================================================================
* OPTIONAL: MERGE EXTERNAL DATA
* ==============================================================================

* If you have deposits data, merge it here:
* merge m:1 protein_id using "deposits_dta.dta", keep(master match)
* replace num_deposits = 0 if missing(num_deposits)
* gen no_deposits = (num_deposits == 0)

* If you have monthly deposits:
* merge 1:1 gene_id ym using "deposits_monthly.dta", keep(master match)

* ==============================================================================
* CREATE GROUP INTERACTIONS FOR HETEROGENEITY
* ==============================================================================

* Create gene group × quarter fixed effects (if gene_group exists)
capture confirm variable gene_group
if !_rc {
    * Clean gene group variable (handle missing and multiple groups)
    gen gene_group_clean = gene_group
    replace gene_group_clean = "missing" if missing(gene_group) | gene_group == ""
    
    * For multiple groups (e.g., "442|454"), take first group
    * You might want to handle this differently based on your research question
    replace gene_group_clean = regexs(1) if regexm(gene_group_clean, "^([0-9]+)")
    
    * Create group × time fixed effects
    egen genegroup_quarter = group(gene_group_clean quarter_num)
}

* ==============================================================================
* QUARTERLY COLLAPSE (if desired)
* ==============================================================================

* If you want quarterly analysis, collapse to protein/gene × quarter
* preserve
* 
* collapse (sum) ///
*     n_papers ///
*     n_top25_y n_top10_y n_top05_y ///
*     n_top25_q n_top10_q n_top05_q ///
*     n_top25_b2 n_top10_b2 n_top05_b2 ///
*     n_newcomer_papers n_veteran_papers ///
*     unique_mesh_count new_mesh_count ///
*     (first) gene_id gene_name protein_existence average_plddt ///
*            high_confidence gene_group_clean, ///
*     by(protein_id quarter_num)
* 
* * Recreate treatment variables after collapse
* gen treat_launch = high_confidence * (quarter_num >= tq(2021q3))
* gen treat_expansion = high_confidence * (quarter_num >= tq(2022q3))

* ==============================================================================
* MAIN ANALYSIS
* ==============================================================================

* Define outcome variables
local outcomes n_papers ///
               n_top25_y n_top10_y n_top05_y ///
               n_top25_q n_top10_q n_top05_q ///
               n_newcomer_papers n_veteran_papers ///
               unique_mesh_count new_mesh_count

* Clear previous estimates
est clear

* Main specification: AlphaFold launch effect
foreach outcome of local outcomes {
    display "Running regression for `outcome'"
    
    * Basic specification with gene and time fixed effects
    qui reghdfe `outcome' treat_launch, ///
        absorb(gene_id ym_seq) ///
        vce(cluster gene_id)
    eststo `outcome'_basic
    
    * With gene group interactions (if available)
    capture confirm variable genegroup_quarter
    if !_rc {
        qui reghdfe `outcome' treat_launch, ///
            absorb(gene_id ym_seq genegroup_quarter) ///
            vce(cluster gene_id)
        eststo `outcome'_groups
    }
}

* ==============================================================================
* RESULTS TABLES
* ==============================================================================

* Basic results table
esttab *_basic, ///
    keep(treat_launch _cons) ///
    label se b(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    compress nogaps ///
    title("Effect of AlphaFold Launch on Research Outcomes") ///
    note("Treatment: High confidence proteins × Post-AlphaFold launch (2021q3+)")

* Results with gene group controls (if available)
capture confirm variable genegroup_quarter
if !_rc {
    esttab *_groups, ///
        keep(treat_launch _cons) ///
        label se b(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
        compress nogaps ///
        title("Effect of AlphaFold Launch (with Gene Group Controls)") ///
        note("Treatment: High confidence proteins × Post-AlphaFold launch (2021q3+)")
}

* ==============================================================================
* ALTERNATIVE SPECIFICATIONS
* ==============================================================================

* Event study specification (quarters around AlphaFold launch)
gen quarters_since_af = quarter_num - tq(2021q3)

* Create event study dummies
forvalues q = -8/8 {
    gen event_`q' = (quarters_since_af == `q') * high_confidence
    if `q' < 0 {
        local pre_dummies "`pre_dummies' event_`q'"
    }
    else if `q' > 0 {
        local post_dummies "`post_dummies' event_`q'"
    }
}

* Event study regression (example with n_papers)
reghdfe n_papers `pre_dummies' `post_dummies', ///
    absorb(gene_id ym_seq) ///
    vce(cluster gene_id)

* ==============================================================================
* HETEROGENEITY ANALYSIS
* ==============================================================================

* By gene group (if available)
capture confirm variable gene_group_clean
if !_rc {
    * Create treatment × gene group interactions
    tab gene_group_clean, gen(group_)
    
    foreach var of varlist group_* {
        gen treat_`var' = treat_launch * `var'
    }
    
    * Regression with group interactions
    reghdfe n_papers treat_launch treat_group_*, ///
        absorb(gene_id ym_seq) ///
        vce(cluster gene_id)
}

* By protein existence level
levelsof protein_existence, local(existence_levels)
foreach level of local existence_levels {
    gen existence_`level' = (protein_existence == "`level'")
    gen treat_existence_`level' = treat_launch * existence_`level'
}

* ==============================================================================
* SUMMARY STATISTICS
* ==============================================================================

* Panel balance check
tab year month, missing

* Treatment balance
bysort high_confidence: summarize average_plddt gene_id

* Outcome summary by treatment status
bysort high_confidence after_af_launch: summarize n_papers n_top25_y unique_mesh_count

display "Analysis complete!"
display "Dataset: final_panel_COMPLETE.dta"
display "Observations: " _N
display "Genes: " `=scalar(_N)/48'
display "Time periods: 2020-2023 (monthly)"
display "Treatment: High AlphaFold confidence (pLDDT > 70)"