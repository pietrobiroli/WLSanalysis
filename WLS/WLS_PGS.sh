# TO RUN THIS FILE:
# $ nohup bash WLS_PGS.sh > WLS_PGS.log &

# This files creates the EA, cog, and noncog PGS for the WLS dataset using the LDpred-inf weights constructed by Travis
# Author: Pietro Biroli
# Date: Feb 2019

######### SECTIONS
# (1) Import the imputed genetic data and clean it (plink)
# (2) Create the polygenic scores using the sum-stats constructed by Travis Mallard


INPUT_GENE='/mnt/department/econ/biroli/geighei/data/WLS/CleanData/WLS_Genotypic_Data'
INPUT_PATH='/mnt/data/Research/cogNoncogGSEM/WLSanalysis'
INPUT_SUMSTAT='/mnt/data/Research/cogNoncogGSEM/3_Summary_statistics/LDpred_sumstats/excl_WLS'

OUTPUT_PATH='/mnt/data/Research/cogNoncogGSEM/PGS/WLS'

cd $OUTPUT_PATH

# printf "\n#-----------------------------------------------------------\n"
# printf "#----- (0) Get the raw data from the source folder and do some QC ----#\n\n"
# # composite quality + informativeness + effective sample size filters (seggested by WLS readme)
# plink1.9 --bfile $INPUT_GENE/Herd_WLS_TOP_subject_level --make-bed --extract $INPUT_GENE/AA_Readme/SNP_composite_filter_effN_extract.txt --out $OUTPUT_PATH/cleanWSL
#
# # all filters + only unrelated individuals
# plink1.9 --bfile $INPUT_GENE/Herd_WLS_TOP_subject_level --make-bed --extract $INPUT_GENE/AA_Readme/SNP_composite_filter_effN_extract.txt --keep $INPUT_GENE/AA_Readme/unrel.keep.txt --out $OUTPUT_PATH/cleanWSL_unrelated


printf "\n#-----------------------------------------------------------\n"
printf "#----- (1) Create the Polygenice gene scores           ----#\n\n"
# column order: rsID allele-code beta-effect-size
head $INPUT_SUMSTAT/COG_excl_WLS/COG_excl_WLS.LDpred-r400_LDpred-inf.txt
# chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta
# chrom_1    752566    rs3094315    G    A    -8.7262e-03    1.5120e-04

## COG
plink1.9 --bfile $INPUT_GENE/cleanWSL_unrelated --score $INPUT_SUMSTAT/COG_excl_WLS/COG_excl_WLS.LDpred-r400_LDpred-inf.txt 3 5 7 center --out COG_PGS_WLS --allow-no-sex

## NONCOG
plink1.9 --bfile $INPUT_GENE/cleanWSL_unrelated --score $INPUT_SUMSTAT/NONCOG_excl_WLS/NONCOG_excl_WLS.LDpred-r400_LDpred-inf.txt 3 5 7 center --out NONCOG_PGS_WLS --allow-no-sex

## EA
plink1.9 --bfile $INPUT_GENE/cleanWSL_unrelated --score $INPUT_SUMSTAT/EA3_excl23andMe_WLS/EA3_excl23andMe_WLS.LDpred-r400_LDpred-inf.txt 3 5 7 center --out EA3_PGS_WLS --allow-no-sex
