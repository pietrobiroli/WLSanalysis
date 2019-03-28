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

*------------ Working directories
if "`c(username)'" == "pbirol" {
   global DIRPGS  "/mnt/data/Research/cogNoncogGSEM/PGS/WLS"
   global DIRDATA "/mnt/department/econ/biroli/geighei/data/WLS/"
   global DIRCODE "/mnt/data/Research/cogNoncogGSEM/analysis/WLS"
}

if "`c(username)'" == "michellerosenberger" {
   global DIRPGS  "/Volumes/g_econ_department$/econ/biroli/geighei/data/WLS/CleanData/PGS"
   global DIRDATA "/Volumes/g_econ_department$/econ/biroli/geighei/data/WLS/"
   global DIRCODE "~/Development/GEIGHEI/WLSanalysis/WLS"
}

cd ${DIRCODE}/output

log using cogNoncog_WLS.log, replace

*------------ Global variables
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

rename fid subject_id   // fid already numeric
rename iid iid_n        // iid already numeric

/*  fid and iid are the same, expect for those that a NA for one of them. I drop them and end up with the sampe of 9012 individuals that is suggested by the WLS

"After the GAC’s standardized QC procedures, including resolution of sample quality and identity, 
genotypes are available on dbGaP for 9,012 unique WLS study participants. 
(Note for 15 pairs of monozygotic twins, only one member of each pair is retained in the unique set of 9,012 participants.)"


gen diff = subject_id - iid_n
sum diff
*/

drop if missing(subject_id)==1
drop if missing(iid_n)==1
* drop fid iid iid_n
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
keep idpriv-rlifewv srbmi-rbmc6  *iq*  z_rh00* z_edeqyr edat16 hidg64 z_ha103rea  ///
     z_jh001rec-z_jg361rer ///personality
     z_jzg01rer hsr* ri001re si001re z_ie008re z_hi101re

destring subject_id , replace
drop if missing(subject_id)==1
duplicates report subject_id if subject_id <.
duplicates drop subject_id, force
merge 1:1 subject_id using "${DIRPGS}/cogNoncogPGS.dta"
drop _merge

*------------ Clean and rename some vars
mvdecode sibcount-rlifewv z_spwiiq_f-z_rh009rec z_ha103rea z_jh001rec-z_jg361rer z_jh001rei-z_jh009rei z_ie008re z_jzg01rer z_hi101re,  mv(-1 = .a \ -2 = .b \ -3 = .r \ -4 = .d \ -5 = .e \ -27 = .f \ -29 = .g)

//Gender
recode z_sexrsp (2 = 0)
rename z_sexrsp male
tab male

//White
recode z_ie008re (2 = 0)
rename z_ie008re white
tab white

//Birth year
rename z_brdxdy yob
replace yob = 1900+yob
replace yob = .w if yob<1900

//Education
rename z_edeqyr eduyears

//Personality
//Note: this variables have more observations than the "phone questionnaire"
rename z_rh001rec extra
rename z_rh003rec openn
rename z_rh005rec neuro
rename z_rh007rec consc
rename z_rh009rec agree
*/

rename z_jn001rei autonomy
rename z_jn010rei envmastery
rename z_jn019rei persgrowth
rename z_jn028rei posrel
rename z_jn037rei purposelife
rename z_jn046rei selfaccept


/* phone questionnaire
rename z_jh001rei extra
rename z_jh032rei openn
rename z_jh025rei neuro
rename z_jh017rei consc
rename z_jh009rei agree
*/


//IQ
rename gwiiq_bm iq

//Attractiveness
rename z_ha103rea attractive

//GPA
rename z_jzg01rer hs_classrank
rename hsrankq    gpa_hsrank        // High school grades percentile rank
rename hsrscorq   gpa_hsranknorm    // High school grades percentile rank-normalized

//Cognition
rename z_hi101re  cog_sim           // WAIS questions in the Cognition-Similarities Module
replace cog_sim = . if cog_sim == -30

rename ri001re    cog_grad          // Total cognition score GRAD (8-items asked)
rename si001re    cog_sib           // Total cognition score SIB (10-items asked) 


/* Cognition variables
z_gi210rec	R5 Letter Fluency: Total number of scored words named.
z_gi310rec	R5 Category Fluency: Total number of scored words named.
z_gi101re	R5 6-item score for the Cognition Similarities section.
z_gi106re	R5 9-item score for the Cognition Similarities section.
z_gi413re	R5 How many correct words did participant remember in immediate recall task?
z_ai101re	R5 6-item score for the Cognition Similarities section.
z_ai413re	R5 How many correct words did R remember in immediate recall task?
z_gi503re	R5 Score for digits ordering task.
z_hi210rec	R6 Letter Fluency: Total number of scored words named.
z_hi310rec	R6 Category Fluency: Total number of scored words named.
z_hi101re	R6 6-item score summarizing Participant's performance on Weschler Adult
            Intelligence Scale (WAIS) questions selected for inclusion in the Cognition - Similarities Module.
z_hi503re	R6 Score for the Digit Ordering Task if Participant never refused an item.
z_hinsraschscore	R6 NS Rasch Scoring of all answered questions
z_hisimplescore	R6 Simple Score for McArdle Task
*/


/* Reverse Raw PGSs beacuse of LDpred
foreach var of varlist *PGS{
	replace `var' = `var'*-1
	}
*/

//Standardize
foreach var of varlist *PGS{
	egen `var'_std = std(`var')
	replace `var' = `var'_std
	drop `var'_std
	}


save "${DIRDATA}/CleanData/WLS_cognoncogPred.dta", replace

} // end CREATEDATA



if $RUNREG==1{
use "${DIRDATA}/CleanData/WLS_cognoncogPred.dta", clear

****************************************************************************
*****************  MULTIVARIATE ANALYSIS - WLS *****************************
****************************************************************************
// Only focuse on race/Origin white
keep if white == 1

gen yob2 = yob^2

// standardize IQ
egen temp = std(iq)
replace iq = temp
drop temp

// Set phenotypes
global PHENO "eduyears iq extra openn neuro consc agree attractive gpa_hsrank cog_sim"


*--------------PGS Correlations (as in CogNonCog_181121.do)
foreach pgs in NONCOG_PGS COG_PGS {
   reg `pgs' ev1-ev20
   predict r`pgs', r 
   egen zr`pgs'=std(r`pgs')
}

corr zrNONCOG_PGS zrCOG_PGS
matrix X = r(C)
putexcel set "${DIRCODE}/output/CogNonCog_WLS.xlsx", sheet(WLSPGSCorr) modify
putexcel a1=matrix(X), names


*--------------Univariate Analysis (as in CogNonCog_181121.do)
// Univariate Program
capture program drop FX 
program define FX
   args y 
   matrix Fx = J(1,6,999)
   capture drop Y
   egen Y = std(`y')
   foreach x in NONCOG_PGS COG_PGS {
      capture drop Z
      egen Z = std(`x')
      reg Y ev1-ev20 male, robust
      local R2 = e(r2)
      reg Y Z ev1-ev20 male, robust 
      matrix A = e(N), _b[Z], _se[Z], e(r2)-`R2', 2*ttail(e(df_r), abs(_b[Z]/_se[Z])), sqrt(e(r2)-`R2')
      matrix rownames A = `x'
      matrix Fx = Fx \ A
      drop Z
      }
   matrix Fx = Fx[2...,1...]
   matrix colnames Fx = N Beta SE DeltaR2 Pval SqrtDeltaR2
   matrix `y' = Fx
   drop Y
end

// Apply Program
foreach pheno in $PHENO {
   di "Univariate Analysis for `pheno'"
   FX `pheno'
}

matrix Fx1 = J(1,12,999)
matrix colnames Fx1 = N Beta_NC SE DeltaR2 Pval SqrtDeltaR2 N Beta_Cog SE DeltaR2 Pval SqrtDeltaR2

foreach pheno in $PHENO {
	matrix `pheno' = `pheno'[1,1...], `pheno'[2,1...] 
	matrix rownames `pheno' = `pheno'
	matrix Fx1 = Fx1 \ `pheno'
	}	
matrix Fx1 = Fx1[2...,1...]	
matrix list Fx1

* matrix list Fx1
putexcel set "${DIRCODE}/output/CogNonCog_WLS.xlsx", sheet(WLSUnivar) modify
putexcel a1 = matrix(Fx1), names


*--------------Multivariate Analysis (as in CogNonCog_181121.do)
// Multivariate Program
capture program drop FXmv 
program define FXmv
   args y 
   preserve
   capture drop Y
   egen Y = std(`y')
   foreach x in NONCOG_PGS COG_PGS {
      egen z`x' = std(`x')
      }
   reg Y zCOG_PGS ev1-ev20 male, robust
   local R2cog = e(r2)
   reg Y zNONCOG_PGS ev1-ev20 male, robust
   local R2nc = e(r2)	
   reg Y zNONCOG_PGS zCOG_PGS ev1-ev20 male, robust 
   matrix A = e(N), _b[zNONCOG_PGS], _se[zNONCOG_PGS], e(r2)-`R2cog', 2*ttail(e(df_r), abs(_b[zNONCOG_PGS]/_se[zNONCOG_PGS])) , sqrt(e(r2)-`R2cog'), ///
               _b[zCOG_PGS], _se[zCOG_PGS], e(r2)-`R2nc', 2*ttail(e(df_r), abs(_b[zCOG_PGS]/_se[zCOG_PGS])), sqrt(e(r2)-`R2nc') 
   matrix rownames A = `y'
   matrix colnames A = N Beta_nc SE_nc DeltaR2_nc Pval_nc SqrtDeltaR2_nc Beta_c SE_c DeltaR2_c Pval_c SqrtDeltaR2_c
   matrix `y' = A
   drop Y
   restore
end

// Apply Program
matrix mvFX = J(1,11,999)

foreach pheno in $PHENO {
	FXmv `pheno'
	matrix mvFX = mvFX \ `pheno'
}

matrix mvFX = mvFX[2...,1...]
matrix colnames mvFX = N Beta_nc SE_nc DeltaR2_nc Pval_nc SqrtDeltaR2_nc Beta_c SE_c DeltaR2_c Pval_c SqrtDeltaR2_c
matrix list mvFX

putexcel set "${DIRCODE}/output/CogNonCog_WLS.xlsx", sheet(WLSMultivar) modify
putexcel a1 = matrix(mvFX), names




*--------------Regressions
foreach yvar in eduyears iq extra openn neuro consc agree attractive {
	des `yvar'
	sum `yvar'
	
	reg `yvar' COG_PGS            yob yob2 male ev1-ev20, cluster(familypriv)
	estimate store cog_`yvar'
	
	reg `yvar'         NONCOG_PGS yob yob2 male ev1-ev20, cluster(familypriv)
	estimate store noncog_`yvar'

	reg `yvar' COG_PGS NONCOG_PGS yob yob2 male ev1-ev20, cluster(familypriv)
	estimate store both_`yvar'
} //end forloop yvar



*--------------Coefficient plot output
*-- Cognitive PGS
coefplot  (both_eduyears , aseq(Years of edu)) ///
          (both_iq , aseq(IQ)) ///
          (both_extra, aseq(Extraversion)) ///
          (both_openn, aseq(Openness)) ///
          (both_neuro, aseq(Neuroticism)) ///
          (both_consc, aseq(Conscientiousness)) ///
          (both_agree, aseq(Agreeableness)) ///
          (both_attractive, aseq(Attractiveness)) ///
   , keep(COG_PGS) ///
   msymbol(D) xline(0) byopts(xrescale) levels(95 90) ciopts(recast(. rcap)) legend(order(1 "95% c.i." 2 "90% c.i." )) graphregion(color(white)) bgcolor(white) ///ciopts(lwidth(2 ..) lcolor(*.2 *.4 ; )) xlabel(-.12(.04).12) 
   aseq swapnames 

graph save   coefplot_cog, replace
graph export coefplot_cog.png, replace

*-- Noncognitive PGS
coefplot  (both_eduyears , aseq(Years of edu)) ///
          (both_iq , aseq(IQ)) ///
          (both_extra, aseq(Extraversion)) ///
          (both_openn, aseq(Openness)) ///
          (both_neuro, aseq(Neuroticism)) ///
          (both_consc, aseq(Conscientiousness)) ///
          (both_agree, aseq(Agreeableness)) ///
          (both_attractive, aseq(Attractiveness)) ///
   , keep(NONCOG_PGS) ///
   msymbol(D) xline(0) byopts(xrescale) levels(95 90) ciopts(recast(. rcap)) legend(order(1 "95% c.i." 2 "90% c.i." )) graphregion(color(white)) bgcolor(white) ///ciopts(lwidth(2 ..) lcolor(*.2 *.4 ; )) xlabel(-.12(.04).12) 
   aseq swapnames 

graph save   coefplot_noncog, replace
graph export coefplot_noncog.png, replace


*-- Both PGS
coefplot  (both_eduyears , aseq(Years of edu)) ///
          (both_iq , aseq(IQ)) ///
          (both_extra, aseq(Extraversion)) ///
          (both_openn, aseq(Openness)) ///
          (both_neuro, aseq(Neuroticism)) ///
          (both_consc, aseq(Conscientiousness)) ///
          (both_agree, aseq(Agreeableness)) ///
          (both_attractive, aseq(Attractiveness)) ///
   , keep(COG_PGS NONCOG_PGS) ///
   msymbol(D) xline(0) byopts(xrescale) levels(95 90) ciopts(recast(. rcap)) legend(order(1 "95% c.i." 2 "90% c.i." )) graphregion(color(white)) bgcolor(white) ///ciopts(lwidth(2 ..) lcolor(*.2 *.4 ; )) xlabel(-.12(.04).12) 
   aseq swapnames 

graph save   coefplot_both, replace
graph export coefplot_both.png, replace


} // end RUNREG

log close
