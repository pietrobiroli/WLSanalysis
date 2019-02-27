* Project:		cogNoncog
* Content:		prediction Analysis
* Author:		Pietro 
* Adapted by:		Dan Belsky
* Date:			February 2019

/* 
This code runs the prediction analysis of the cog and noncog PGS in the WLS data
*/

********************************************************************************
*********************************** PREAMBLE ***********************************
********************************************************************************
capture log close
clear all
est clear
set more off
set emptycells drop
set matsize 10000
set maxvar 20000

global DIRPGS  "/mnt/data/Research/cogNoncogGSEM/PGS/WLS"
global DIRDATA "/mnt/department/econ/biroli/geighei/data/WLS/"
global DIRCODE "/mnt/data/Research/cogNoncogGSEM/analysis/WLS"


cd ${DIRCODE}/output

log using cogNoncog_WLS.log, replace

global CREATEDATA = 1
global RUNREG     = 1



if $CREATEDATA==1{
********************************************************************************
*******************************MERGE THE DATA **********************************
********************************************************************************
*------------ Import the PGS created in plink
import delimited "${DIRPGS}/COG_PGS_WLS.profile", clear delimiter(space, collapse)
drop v1 pheno cnt cnt2 
rename score COG_PGS
compress
save "${DIRPGS}/cogNoncogPGS.dta", replace

import delimited "${DIRPGS}/NONCOG_PGS_WLS.profile", clear delimiter(space, collapse)
drop v1 pheno cnt cnt2 
rename score NONCOG_PGS
compress
merge 1:1 fid iid using "${DIRPGS}/cogNoncogPGS.dta"
drop _merge
save "${DIRPGS}/cogNoncogPGS.dta", replace

import delimited "${DIRPGS}/EA3_PGS_WLS.profile", clear delimiter(space, collapse)
drop v1 pheno cnt cnt2 
rename score EA3_PGS
compress
merge 1:1 fid iid using "${DIRPGS}/cogNoncogPGS.dta"
drop _merge

duplicates report fid iid
destring fid, force gen(subject_id)
destring iid, force gen(iid_n)

/*  fid and iid are the same, expect for those that a NA for one of them. I drop them and end up with the sampe of 9012 individuals that is suggested by the WLS

"After the GAC’s standardized QC procedures, including resolution of sample quality and identity, 
genotypes are available on dbGaP for 9,012 unique WLS study participants. 
(Note for 15 pairs of monozygotic twins, only one member of each pair is retained in the unique set of 9,012 participants.)"


gen diff = subject_id - iid_n
sum diff
*/

drop if missing(subject_id)==1
drop if missing(iid_n)==1
drop fid iid iid_n
save "${DIRPGS}/cogNoncogPGS.dta", replace

*------------ Merge with the principal components
import delimited "${DIRDATA}/RawData/WLS_Genotypic_Data/AA_Readme/Principal_components.csv", clear delimiter(comma, collapse)
merge 1:1 subject_id using "${DIRPGS}/cogNoncogPGS.dta"
drop _merge
save "${DIRPGS}/cogNoncogPGS.dta", replace

*------------ Merge with the information data
import delimited "${DIRDATA}/RawData/WLS_Genotypic_Data/AA_Readme/Sample_analysis.csv", clear delimiter(comma, collapse)
destring subject_id, replace force 
drop if missing(subject_id)==1
drop sample_id 
merge 1:1 subject_id using "${DIRPGS}/cogNoncogPGS.dta"
drop _merge
save "${DIRPGS}/cogNoncogPGS.dta", replace


*------------ Merge with the phenotypic data
use "${DIRDATA}/RawData/WLS_Survey_PhenotypicLong-formData-13_06/wls_plg_13_06.dta", clear
keep idpriv-rlifewv srbmi-rbmc6  *iq*  z_rh00*

destring subject_id , replace
drop if missing(subject_id)==1
duplicates report subject_id if subject_id <.
duplicates drop subject_id, force
merge 1:1 subject_id using "${DIRPGS}/cogNoncogPGS.dta"
drop _merge


*------------ Clean and rename some vars
recode z_sexrsp (2 = 0)
rename z_sexrsp male
tab male

rename z_brdxdy yob
replace yob = 1900+yob
replace yob = .w if yob<1900


//Personality
rename z_rh001rec extra
rename z_rh003rec openn
rename z_rh005rec neuro
rename z_rh007rec consc
rename z_rh009rec agree

//IQ
rename gwiiq_bm iq
	
* Reverse Raw PGSs beacuse of LDpred
foreach var of varlist *PGS{
	replace `var' = `var'*-1
	}
*

//Standardize
foreach var of varlist *PGS{
	egen `var'_std = std(`var')
	replace `var' = `var'_std
	drop `var'_std
	}

mvdecode iq extra openn neuro consc agree, mv(-3 = .r)

save "${DIRDATA}/CleanData/WLS_cognoncogPred.dta", replace
} // end CREATEDATA



if $RUNREG==1{
use "${DIRDATA}/CleanData/WLS_cognoncogPred.dta", clear

****************************************************************************
*****************  MULTIVARIATE ANALYSIS - WLS *****************************
****************************************************************************
gen yob2 = yob^2

foreach yvar in iq extra openn neuro consc agree{
	des `yvar'
	sum `yvar'
	reg `yvar' COG_PGS NONCOG_PGS yob yob2 male ev1-ev20, cluster(familypriv)
}

} // end RUNREG

log close
