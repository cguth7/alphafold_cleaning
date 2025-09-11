cd "/Users/maxguthmann/Dropbox/Matteo_Danilo_AlphaFold2"



************************************************************
************ PUBTATOR: keep disease-revelant  pubs *********
************************************************************


// 1. Load the pubtator data
import delimited "rawdata/disease2pubtator3.txt", delimiter(tab) clear
rename (v1 v2 v3 v4 v5) (PMID Type Mesh_ID Disease_Name Source)

// Keep only PMID and Mesh_ID
keep PMID Mesh_ID

// 2. Flag OMIM rows
gen is_omim = strpos(Mesh_ID, "OMIM")>0

// 3. Strip off "MESH:" if present
replace Mesh_ID = subinstr(Mesh_ID, "OMIM:", "", .)
replace Mesh_ID = subinstr(Mesh_ID, "MESH:", "", .)

// 4. Save for merging
save "filedata/disease_mesh.dta", replace


//////////////////////// culll out the non 3/4/5 levels


use "rawdata/charlie/mesh_depth_3_4_5_no_explode.dta", clear
rename DescriptorUI Mesh_ID   // adjust if your column is named differently
keep Mesh_ID
duplicates drop
sort Mesh_ID
tempfile meshlist
save `meshlist'

// 6. Re-load pubtator, merge with meshlist
use "filedata/disease_mesh.dta", clear
merge m:1 Mesh_ID using `meshlist'

// 7. Keep either: (a) OMIM rows, or (b) MeSH that matched depth list
keep if is_omim==1 | _merge==3
drop _merge

save "filedata/pmid_mesh_depth345_or_omim.dta", replace


// 8. Drop to unique PMIDs
keep PMID
duplicates drop PMID, force
sort PMID

// 9. Save final set
save "filedata/unique_pmids_with_disease_depth345_or_omim.dta", replace

************************************************************
************ PUBTATOR: keep gene-relevant pubs *************
************ among the disease-revelant  pubs **************
************************************************************


import delimited "rawdata/gene2pubtator3.txt", delimiter(tab) clear

// rename the columns
rename (v1 v2 v3 v4 v5) (PMID Type Gene_ID Gene_Name Source) 

// probably don't need to
sort PMID

// inner merge, we have unique ids in the using file not in the file currently read into stata
merge m:1 PMID using "filedata/unique_pmids_with_disease_depth345_or_omim.dta"
keep if _merge == 3
drop _merge

save "filedata/gene_in_disease_pubs_3_4_5.dta", replace




//////////////////
//// START HERE IF YOU DON'T WANT IT TO TAKE FOREVER /////
//////////////////

************************************************************
************ Master: MERGE IN MASTER TO GET REDUCED GENE *************
************ among the disease-revelant  pubs **************
************************************************************

///keeps freezing here so i struggle.


* unique genes from master
use "filedata/protein_data_master_dta.dta", clear
rename geneid Gene_ID
tostring Gene_ID, replace force
keep Gene_ID
duplicates drop
tempfile keepgenes
save `keepgenes'

* filter paper–gene to those genes
use "/Users/maxguthmann/Downloads/Development/Work/unique_diseases/gene_in_disease_cleaned.dta", clear
capture confirm string variable Gene_ID
if _rc tostring Gene_ID, replace force

merge m:1 Gene_ID using `keepgenes'
keep if _merge==3
drop _merge

// save "filedata/gene_in_disease_pubs.dta", replace


************************************************************
************ DATES from NIH DATA  *************
************************************************************

use "filedata/gene_in_disease_pubs.dta", clear


// need to rename PMID to pmid for merge
rename PMID pmid

// add the dates to the database -> we lose ~14m out of 63m. 
merge m:1 pmid using "rawdata/complete_pmid_dates.dta"
keep if _merge==3
drop _merge

// only keep papers between 2020 and 2023. 

//check if anything is wrong although assumption is no null dates from the data set. (confirmed but left this in)
assert inrange(month,1,12) if !missing(month)
drop if missing(year) | missing(month)


//create cat vars for month date in statas datetime (27,234,856 observations deleted) (also 2024 less than 1% of the density)
gen ym = ym(year, month)
format ym %tm
keep if inrange(ym, ym(2020,1), ym(2023,12))

save "filedata/gene_in_disease_pubs_2020_2023_MESH_3_4_5.dta", replace


************************************************************
************ ICITE: get bibliographic details  *************
************************************************************


/// create percentiles from bibliographic data

use "filedata/protein_data_master_dta.dta", clear
rename geneid Gene_ID
// tostring Gene_ID, replace force
keep Gene_ID
duplicates drop
tempfile keepgenes
save `keepgenes'


/////////////NOTE THIS IS A DISCONTINOUS BREAK WITH A PYTHON DATA INSERT FIX
use "filedata/gene_in_disease_pubs_2020_2023_MESH_all_levels_3_4_5.dta", clear
drop if strpos(Gene_ID, ";") > 0
destring Gene_ID, replace


merge m:1 Gene_ID using `keepgenes'
keep if _merge==3
drop _merge

save "filedata/gene_in_disease_pubs_filtered_genes_2020_2023_MESH_all_levels_3_4_5.dta", replace



/// get a temp file and use it 3 times

* 0) Start: merge iCite
use "filedata/gene_in_disease_pubs_filtered_genes_2020_2023_MESH_all_levels_3_4_5.dta", clear
merge m:1 pmid using "rawdata/charlie/icite_cites_join.dta"
keep if _merge==1 | _merge==3
drop _merge
replace icite_cites = 0 if missing(icite_cites)

* 1) Paper-level base: KEEP year month ym from master
tempfile papers yr qrtr b2 final_paper_flags
preserve
    collapse (firstnm) year month ym icite_cites, by(pmid)
    * quarter & bi-month derived from existing ym/month
    gen quarter_num = qofd(dofm(ym))
    format quarter_num %tq

    gen b2_bin      = floor((month-1)/2)         // 0..5 within year
    gen b2_index    = year*6 + b2_bin
    gen b2_start_ym = ym(year, 2*b2_bin + 1)
    format b2_start_ym %tm

    save `papers'
restore

///////





// gen quarter and bi-monthly

*******************************************************
*** 2) YEAR bins → hi25 hi10 hi05  (pctile cutoffs)
*******************************************************
use `papers', clear
bys year: egen p75 = pctile(icite_cites), p(75)
bys year: egen p90 = pctile(icite_cites), p(90)
bys year: egen p95 = pctile(icite_cites), p(95)

gen hi25 = icite_cites >= p75
gen hi10 = icite_cites >= p90
gen hi05 = icite_cites >= p95
label var hi25 "Top 25% cites (within year)"
label var hi10 "Top 10% cites (within year)"
label var hi05 "Top 5% cites (within year)"

keep pmid hi25 hi10 hi05
save `yr'

*******************************************************
*** 3) QUARTER bins → hi25_q hi10_q hi05_q
*******************************************************
use `papers', clear
bys quarter_num: egen p75_q = pctile(icite_cites), p(75)
bys quarter_num: egen p90_q = pctile(icite_cites), p(90)
bys quarter_num: egen p95_q = pctile(icite_cites), p(95)

gen hi25_q = icite_cites >= p75_q
gen hi10_q = icite_cites >= p90_q
gen hi05_q = icite_cites >= p95_q
label var hi25_q "Top 25% (within quarter)"
label var hi10_q "Top 10% (within quarter)"
label var hi05_q "Top 5%  (within quarter)"

keep pmid hi25_q hi10_q hi05_q
save `qrtr'

*******************************************************
*** 4) BI-MONTH bins → hi25_b2 hi10_b2 hi05_b2
*******************************************************
use `papers', clear
bys b2_index: egen p75_b2 = pctile(icite_cites), p(75)
bys b2_index: egen p90_b2 = pctile(icite_cites), p(90)
bys b2_index: egen p95_b2 = pctile(icite_cites), p(95)

gen hi25_b2 = icite_cites >= p75_b2
gen hi10_b2 = icite_cites >= p90_b2
gen hi05_b2 = icite_cites >= p95_b2
label var hi25_b2 "Top 25% (within bi-month)"
label var hi10_b2 "Top 10% (within bi-month)"
label var hi05_b2 "Top 5%  (within bi-month)"

keep pmid hi25_b2 hi10_b2 hi05_b2
save `b2'

*******************************************************
*** 5) Smash flags together + keep meta cols
*******************************************************
use `papers', clear
keep pmid year month ym quarter_num b2_index b2_start_ym icite_cites
merge 1:1 pmid using `yr',   nogen
merge 1:1 pmid using `qrtr', nogen
merge 1:1 pmid using `b2',   nogen
save "filedata/charlie/icite_highcites_paperlevel_by_year_quarter_bi_monthly.dta", replace
save `final_paper_flags'

*******************************************************
*** 6) Merge back to gene–paper rows and save final
*******************************************************
use "filedata/gene_in_disease_pubs_filtered_genes_2020_2023_MESH_all_levels_3_4_5.dta", clear
merge m:1 pmid using `final_paper_flags'
keep if _merge==3
drop _merge

save "filedata/gene_in_disease_pubs_filtered_genes_2020_2023_MESH_all_levels_3_4_5_with_citation_percentile_by_year_quarterly_bi_monthly.dta", replace



************************************************************
************ Aggregations  *************
************************************************************


///// citation aggregations



use "filedata/gene_in_disease_pubs_filtered_genes_2020_2023_MESH_all_levels_3_4_5_with_citation_percentile_by_year_quarterly_bi_monthly.dta", clear


/// 20k missing at one point?, from pubtator, so small it doesn't matter
drop if missing(Gene_ID)


// 1) One row per Gene_ID × ym × pmid; carry max() of each flag
bysort Gene_ID ym pmid: egen hi25_any      = max(hi25)       // year-based
bysort Gene_ID ym pmid: egen hi10_any      = max(hi10)
bysort Gene_ID ym pmid: egen hi05_any      = max(hi05)

bysort Gene_ID ym pmid: egen hi25_q_any    = max(hi25_q)     // quarter-based
bysort Gene_ID ym pmid: egen hi10_q_any    = max(hi10_q)
bysort Gene_ID ym pmid: egen hi05_q_any    = max(hi05_q)

bysort Gene_ID ym pmid: egen hi25_b2_any   = max(hi25_b2)    // bi-month-based
bysort Gene_ID ym pmid: egen hi10_b2_any   = max(hi10_b2)
bysort Gene_ID ym pmid: egen hi05_b2_any   = max(hi05_b2)

// keep a single row per triplet after computing the *_any flags
bysort Gene_ID ym pmid: keep if _n == 1

// 2) Collapse to Gene_ID × ym with counts for each definition
collapse (count) n_papers = pmid ///
         (sum) ///
             n_top25_y  = hi25_any   n_top10_y  = hi10_any   n_top05_y  = hi05_any   ///
             n_top25_q  = hi25_q_any n_top10_q  = hi10_q_any n_top05_q  = hi05_q_any ///
             n_top25_b2 = hi25_b2_any n_top10_b2 = hi10_b2_any n_top05_b2 = hi05_b2_any, ///
         by(Gene_ID ym)
		 
format ym %tm
sort Gene_ID ym

save "filedata/gene_month_citation_counts_latest_by_year_quarter_bi_monthly.dta", replace


**************************************************************************************
************ Aggregations. ----> loading in unique/ new_unique data from python  *************   /Users/maxguthmann/Downloads/Development/Work/unique_diseases/creating_unique_set.ipynb
**************************************************************************************


use "filedata/gene_month_citation_counts_latest_by_year_quarter_bi_monthly.dta", clear

//// annoying but i guess change dt to merge > reason is because we have (a very very small amount of data that is like10661;10365;51274;9314;688;1316;8609;11279;687;7071;8462;11278;51621;136259;28999;83855;128209;105378  ).  <- we can just drop later if matteo doesn't care, or cut the strings but if we go lower than 50 we get duplciates

// recast str244 Gene_ID, force
// tostring Gene_ID, replace
// tostring ym, replace


//// change to relative, some issues with dropbox syncing. 

merge 1:1 Gene_ID ym using "/Users/maxguthmann/Dropbox/Matteo_Danilo_AlphaFold2/filedata/charlie/gene_month_unique_mesh_2015.dta", keep(master match)
drop _merge

merge 1:1 Gene_ID ym using "/Users/maxguthmann/Dropbox/Matteo_Danilo_AlphaFold2/filedata/charlie/gene_month_new_unique_diseases_2015_no_level_5.dta", keep(master match)
replace new_mesh_count = 0 if missing(new_mesh_count)
drop _merge

///clean up 

drop ym_date
// destring Gene_ID, replace

save "filedata/gene_month_citations_and_unique_disease_latest_by_year_quarter_bi_monthly_2015.dta", replace



// maybe : downsize the annoying gene id string. or we could just intify


************************************************************
************ Intersection Variable *************                    ----> need to do this in python, will get back to it. 
************************************************************






//// drop if semi colon, merge onto main, balance the panel. 

use "filedata/protein_data_master_dta.dta", clear
rename geneid Gene_ID
// tostring Gene_ID, replace force
keep Gene_ID protein_id gene_name protein_existence
duplicates drop Gene_ID, force
tempfile keepgenes
save `keepgenes'

use "rawdata/genes_groups.dta", clear
rename ncbi_entrez_gene_id Gene_ID
keep Gene_ID GenegroupID
duplicates drop Gene_ID, force
tempfile genegroups
save `genegroups'

use "filedata/gene_month_citations_and_unique_disease_latest_by_year_quarter_bi_monthly_2015.dta", clear
merge m:1 Gene_ID using `keepgenes'

keep if _merge==3 | _merge==2
drop _merge


merge m:1 Gene_ID using `genegroups'

keep if _merge==3 | _merge==1
drop _merge

//// make it into a panel

// Ensure one row per gene-month before filling
bys Gene_ID ym: keep if _n==1

// Create the full grid of gene x month
fillin Gene_ID ym
drop if missing(ym)


/// need to merge back in to fill the gene names/ protein id for everthing
merge m:1 Gene_ID using `keepgenes', update nogenerate
// *** FIX ENDS HERE ***

// Replace missing citation counts with 0 for the new rows
// Replace missing counts with 0
foreach v in n_papers                         ///
             n_top25_y n_top10_y n_top05_y    ///
             n_top25_q n_top10_q n_top05_q    ///
             n_top25_b2 n_top10_b2 n_top05_b2 ///
             unique_mesh_count new_mesh_count {
    replace `v' = 0 if missing(`v')
}




///// merge in the gene families






// // Use encode to create a safe numeric ID for panel operations
// encode Gene_ID, gen(gene_id_numeric)
//
// // Verify the panel is balanced
// xtset gene_id_numeric ym
// distinct gene_id_numeric
// local unique_genes = r(ndistinct)
// assert _N == `unique_genes' * 48

// Save the final, complete, and balanced panel
/// probably should rename this
save "filedata/final_balanced_panel_inclusive_by_year_quarter_bi_month_2015_with_groups.dta", replace