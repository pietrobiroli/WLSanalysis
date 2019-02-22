global user db3275
//global user dbelsky
//********************************************************************************//
//********************************************************************************//
//DUNEDIN
//********************************************************************************//
//********************************************************************************//
use "/Users/$user/Box/Belsky/NonCogPGS/PRSiceScores/NonCogLDPScores_Dunedin.dta", clear
merge 1:1 snum using "/Users/$user/Box/DPPP_genomics/Genomics_OUT/PGSFiles/PCs/DunedinPCs100.dta", keepus(zpc*) nogen

merge 1:1 snum using "/Users/$user/Box/Belsky/NonCogPGS/NZ_phenodata/PGS_nonCog_NZ_10Sept.dta", keepus(lscuw311-WFSIQ38) nogen
merge 1:1 snum using "/Users/$user/Box/Belsky/CerebCortex2018/Dunedin Phenotypes/IQ3to38.dta", nogen
merge 1:1 snum using "/Users/$user/Box/Belsky/PNAS2018/Code/SocMob_Dunedin_Archive.dta", nogen keepus(Educ263238 sex)	

gen selfcon = lsc*-1

//Personality
rename BF_ex extr
rename BF_ag agre
rename BF_con CON
rename BF_neu neur
rename BF_opn open

//Childhood IQ
foreach x in inf sim ari voc pic blk obj cod {
	foreach y in 7 9 11 13{ 
		egen Z`y' = std(`x'_ma`y')
		}
	egen `x' = rowmean(Z7 Z9 Z11 Z13)
	drop Z7 Z9 Z11 Z13 
	}
	
//Reverse Raw PGSs
foreach x in pgsEA3LDPINF pgsCPLDPINF pgsNONCOGLDPINF pgsCOGLDPINF{
	replace `x' = `x'*-1
	}

save "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_Dunedin_181121.dta", replace

//******************************************************************************//
// PGS Correlations
//******************************************************************************//
use "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_Dunedin_181121.dta", clear

foreach x in pgsEA3LDPINF pgsCPLDPINF pgsNONCOGLDPINF pgsCOGLDPINF {
	reg `x' zpc*
	predict r`x', r
	egen zr`x'=std(r`x')
	}
corr zrpgs*
matrix X = r(C)
putexcel set "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_181121.xlsx", sheet(DunedinPGSCorr) modify
putexcel a1=matrix(X), names


//****************************************************************************//
//****************************************************************************//
//UNIVARIATE ANALYSIS - DUNEDIN
//****************************************************************************//
//****************************************************************************//
	//UNIVARIATE PROGRAM - DUNEDIN 
capture program drop FX 
program define FX
args y 
matrix Fx = J(1,6,999)
capture drop Y
egen Y = std(`y')
foreach x in pgsEA3LDPINF pgsCPLDPINF pgsNONCOGLDPINF pgsCOGLDPINF {
	capture drop Z
	egen Z = std(`x')
	reg Y zpc* sex, robust
	local R2 = e(r2)
	reg Y Z zpc* sex , robust 
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
//****************************************************************************//
#delimit ;
foreach PH in 
	inf voc sim ari blk pic obj cod
	wfsiq713 wviq713 wpiq713 lscuw311
	inf_ss38 voc_ss38 sim_ss38 ar_ss38 bd_ss38 mr_ss38 pc_ss38 ds_ss38 ss_ss38 dsc_ss38
	WFSIQ38 vci38 pri38 wmi38 psi38
	open CON extr agre neur 
	{ ; #delimit cr
	FX `PH'
	}
	
matrix Fx1 = J(1,24,999)
matrix colnames Fx1 = N Beta_EA SE DeltaR2 Pval SqrtDeltaR2 N Beta_CP SE DeltaR2 Pval SqrtDeltaR2 N Beta_NC SE DeltaR2 Pval SqrtDeltaR2 N Beta_Cog SE DeltaR2 Pval SqrtDeltaR2
#delimit ;
foreach PH in 
	inf voc sim ari blk pic obj cod
	wfsiq713 wviq713 wpiq713 lscuw311
	inf_ss38 voc_ss38 sim_ss38 ar_ss38 bd_ss38 mr_ss38 pc_ss38 ds_ss38 ss_ss38 dsc_ss38
	WFSIQ38 vci38 pri38 wmi38 psi38
	open CON extr agre neur 
	{ ; #delimit cr
	matrix `PH' = `PH'[1,1...], `PH'[2,1...], `PH'[3,1...] , `PH'[4,1...] 
	matrix rownames `PH' = `PH'
	matrix Fx1 = Fx1 \ `PH'
	}	
matrix Fx1 = Fx1[2...,1...]	
matrix list Fx1

putexcel set "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_181121.xlsx", sheet(DunedinUnivar) modify
putexcel a1 = matrix(Fx1), names
//********************************************************************************//


//****************************************************************************//
//****************************************************************************//
//MULTIVARIATE ANALYSIS - DUNEDIN 
//****************************************************************************
//****************************************************************************//
	//MULTIVARIATE PROGRAM - DUNEDIN 
capture program drop FXmv 
program define FXmv
args y 
	preserve
	capture drop Y
	egen Y = std(`y')
	foreach x in pgsNONCOGLDPINF pgsCOGLDPINF{
		egen z`x' = std(`x')
		}
	reg Y zpgsCOGLDPINF zpc* sex, robust
	local R2cog = e(r2)
	reg Y zpgsNONCOGLDPINF zpc* sex, robust
	local R2nc = e(r2)	
	reg Y zpgsNONCOGLDPINF zpgsCOGLDPINF zpc* sex , robust 
	matrix A = e(N), _b[zpgsNONCOGLDPINF], _se[zpgsNONCOGLDPINF], e(r2)-`R2cog', 2*ttail(e(df_r), abs(_b[zpgsNONCOGLDPINF]/_se[zpgsNONCOGLDPINF])) , sqrt(e(r2)-`R2cog'), ///
					_b[zpgsCOGLDPINF], _se[zpgsCOGLDPINF], e(r2)-`R2nc', 2*ttail(e(df_r), abs(_b[zpgsCOGLDPINF]/_se[zpgsCOGLDPINF])), sqrt(e(r2)-`R2nc') 
	matrix rownames A = `y'
	matrix colnames A = N Beta_nc SE_nc DeltaR2_nc Pval_nc SqrtDeltaR2_nc Beta_c SE_c DeltaR2_c Pval_c SqrtDeltaR2_c
	matrix `y' = A
	drop Y
	restore
end
//****************************************************************************//
matrix mvFX = J(1,11,999)
#delimit ;
foreach PH in 
	inf voc sim ari blk pic obj cod
	wfsiq713 wviq713 wpiq713 lscuw311
	inf_ss38 voc_ss38 sim_ss38 ar_ss38 bd_ss38 mr_ss38 pc_ss38 ds_ss38 ss_ss38 dsc_ss38
	WFSIQ38 vci38 pri38 wmi38 psi38
	open CON extr agre neur 
	{ ; #delimit cr
	FXmv `PH'
	matrix mvFX = mvFX \ `PH'
	}
matrix mvFX = mvFX[2...,1...]
matrix colnames mvFX = N Beta_nc SE_nc DeltaR2_nc Pval_nc SqrtDeltaR2_nc Beta_c SE_c DeltaR2_c Pval_c SqrtDeltaR2_c
matrix list mvFX

putexcel set "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_181121.xlsx", sheet(DunedinMultivar) modify
putexcel a1 = matrix(mvFX), names
	
//********************************************************************************//




//********************************************************************************//
//********************************************************************************//
//ERISK
//********************************************************************************//
//********************************************************************************//
use "/Users/$user/Box/DPPP_genomics/Genomics_OUT/PGSFiles/PCs/ERisk/ZygosityUpdate_2018.dta", clear
merge 1:1 atwinid using  "/Users/$user/Box/Belsky/NonCogPGS/PRSiceScores/NonCogLDPScores_ERisk.dta", nogen
merge 1:1 atwinid using "/Users/$user/Box/DPPP_genomics/Genomics_OUT/PGSFiles/PCs/ERiskTwinPCs.dta", keepus(zpc*) nogen
merge 1:1 atwinid using "/Users/$user/Box/Belsky/NonCogPGS/Erisk_phenodata/Dan_20Sep18.dta", nogen
merge 1:1 atwinid using "/Users/$user/Box/Belsky/PNAS2018/Code/SocMob_ERisk_Archive.dta", keepus(educachve) nogen

//Personality
rename bfioe18 open
rename bfice18 CON
rename bfiee18 extr
rename bfiae18 agre
rename bfine18 neur


//Reverse Raw PGSs
foreach x in pgsEA3LDPINF pgsCPLDPINF pgsNONCOGLDPINF pgsCOGLDPINF{
	replace `x' = `x'*-1
	}

save "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_ERisk_181121.dta", replace



//******************************************************************************//
// PGS Correlations
//******************************************************************************//
use "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_ERisk_181121.dta", clear

foreach x in pgsEA3LDPINF pgsCPLDPINF pgsNONCOGLDPINF pgsCOGLDPINF {
	reg `x' zpc*
	predict r`x', r
	egen zr`x'=std(r`x')
	}
corr zrpgs*
matrix X = r(C)
putexcel set "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_181121.xlsx", sheet(ERiskPGSCorr) modify
putexcel a1=matrix(X), names


//****************************************************************************//
//UNIVARIATE ANALYSIS - ERISK 
//****************************************************************************//
	//UNIVARIATE PROGRAM - ERISK 
capture program drop FX 
program define FX
args y 
matrix Fx = J(1,6,999)
capture drop Y
egen Y = std(`y')
foreach x in pgsEA3LDPINF pgsCPLDPINF pgsNONCOGLDPINF pgsCOGLDPINF{
	capture drop Z
	egen Z = std(`x')
	reg Y zpc* sampsex, robust
	local R2 = e(r2)
	reg Y Z zpc* sampsex , robust cluster(familyid)
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
//****************************************************************************//
#delimit ;
foreach PH in 
	iq18e infe_ss18 mre_ss18 dsce_ss18
	rvpapre18 rvpmlte18 rvptface18 swmstae18 swmteae18 swmmlre18 sspsple18 ssprsle
	lowsc510e
	open CON extr agre neur 
	{ ; #delimit cr
	FX `PH'
	}	
matrix Fx1 = J(1,24,999)
matrix colnames Fx1 = N Beta_EA SE DeltaR2 Pval SqrtDeltaR2 N Beta_CP SE DeltaR2 Pval SqrtDeltaR2 N Beta_NC SE DeltaR2 Pval SqrtDeltaR2
#delimit ;
foreach PH in 
	iq18e infe_ss18 mre_ss18 dsce_ss18
	rvpapre18 rvpmlte18 rvptface18 swmstae18 swmteae18 swmmlre18 sspsple18 ssprsle
	lowsc510e
	open CON extr agre neur
	{ ; #delimit cr
	matrix `PH' = `PH'[1,1...], `PH'[2,1...], `PH'[3,1...], `PH'[4,1...]
	matrix rownames `PH' = `PH'
	matrix Fx1 = Fx1 \ `PH'
	}	
matrix Fx1 = Fx1[2...,1...]	
matrix list Fx1
putexcel set "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_181121.xlsx", sheet(ERiskUnivar) modify
putexcel a1 = matrix(Fx1), names

//****************************************************************************//
//MULTIVARIATE ANALYSIS - ERISK 
//****************************************************************************//
	//MULTIVARIATE PROGRAM - ERISK 
capture program drop FXmv 
program define FXmv
args y 
	preserve
	capture drop Y
	egen Y = std(`y')
	foreach x in pgsNONCOGLDPINF pgsCOGLDPINF{
		egen z`x' = std(`x')
		}
	reg Y zpgsCOGLDPINF zpc* sampsex, robust
	local R2cog = e(r2)
	reg Y zpgsNONCOGLDPINF zpc* sampsex, robust
	local R2nc = e(r2)	
	reg Y zpgsNONCOGLDPINF zpgsCOGLDPINF zpc* sampsex , robust 
	matrix A = e(N), _b[zpgsNONCOGLDPINF], _se[zpgsNONCOGLDPINF], e(r2)-`R2cog', 2*ttail(e(df_r), abs(_b[zpgsNONCOGLDPINF]/_se[zpgsNONCOGLDPINF])) , sqrt(e(r2)-`R2cog'), ///
					_b[zpgsCOGLDPINF], _se[zpgsCOGLDPINF], e(r2)-`R2nc', 2*ttail(e(df_r), abs(_b[zpgsCOGLDPINF]/_se[zpgsCOGLDPINF])), sqrt(e(r2)-`R2nc') 
	matrix rownames A = `y'
	matrix colnames A = N Beta_nc SE_nc DeltaR2_nc Pval_nc SqrtDeltaR2_nc Beta_c SE_c DeltaR2_c Pval_c SqrtDeltaR2_c
	matrix `y' = A
	drop Y
	restore
end
//****************************************************************************//
matrix mvFX = J(1,11,999)
#delimit ;
foreach PH in 
	iq18e infe_ss18 mre_ss18 dsce_ss18
	rvpapre18 rvpmlte18 rvptface18 swmstae18 swmteae18 swmmlre18 sspsple18 ssprsle
	lowsc510e
	open CON extr agre neur 
	{ ; #delimit cr
	FXmv `PH'
	matrix mvFX = mvFX \ `PH'
	}
matrix mvFX = mvFX[2...,1...]
matrix colnames mvFX = N Beta_nc SE_nc DeltaR2_nc Pval_nc SqrtDeltaR2_nc Beta_c SE_c DeltaR2_c Pval_c SqrtDeltaR2_c
matrix list mvFX
putexcel set "/Users/$user/Box/Belsky/NonCogPGS/Code/CogNonCog_181121.xlsx", sheet(ERiskMultivar) modify
putexcel a1 = matrix(mvFX), names
//********************************************************************************//




