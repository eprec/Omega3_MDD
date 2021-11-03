#######################################################################################################################
## OMEGA 3 MDD MR SCRIPT
#######################################################################################################################

rm(list=ls(all=TRUE)) #empties your R environment
setwd("/Omega 3")

######################################################################################################################################
#1. Install packages
######################################################################################################################################

{
  #Load packages - you will have to install the packages first time you run the script
  
  #install.packages("devtools")
  library(devtools)
  #install_github("MRCIEU/TwoSampleMR") #re-run this to update MR BASE
  library(TwoSampleMR)
  #install.packages("ggplot2")
  library(ggplot2)
  #install.packages("knitr")
  library(knitr)
  #install_github('qingyuanzhao/mr.raps')
  library(mr.raps)
  #install.packages("markdown")
  library(markdown) 
  #install.packages("MRInstruments")
  library(MRInstruments)
  #install.packages("readxl")
  library("readxl")
  library("dplyr")
  library(tibble)
  library(markdown) 
  library("dplyr")
  library(tibble)
  library(MendelianRandomization)
  library(xlsx)
  #install.packages("haven")
  library(haven)
  #install.packages("xlsx")
  library("tidyverse")
  library(ieugwasr)
  
  install.packages("phewas")
  library(phewas)
  
  devtools::install_github("explodecomputer/genetics.binaRies")
  genetics.binaRies::get_plink_binary()
  
  ao <- available_outcomes()
}





ao <- available_outcomes()
gwas_catalog<- gwas_catalog

#######################################################################################################################
##1.  Extract instrument sets for each exposure from each of the following UKBB GWAS studies from IEU Open GWAS project:
#######################################################################################################################
#1.Total Omega 3 
#2.DHA in UKBB
#3.Omega 3 % 
#4 Total Omega 6
#5 Linoleic Acid

exposures <- c("met-d-Omega_3", "met-d-DHA","met-d-Omega_3_pct", "met-d-Omega_6", "met-d-LA")

for (i in 1:5){
  exposure_dat<-extract_instruments(exposures[i])
  exposure_dat<- format_data(exposure_dat, type = "exposure", header = TRUE,
              phenotype_col = "exposure", snp_col = "SNP", beta_col = "beta.exposure", samplesize_col = "samplesize.exposure",
              se_col = "se.exposure", eaf_col = "eaf.exposure", effect_allele_col = "effect_allele.exposure",
              other_allele_col = "other_allele.exposure", pval_col = "pval.exposure")
  exposure_dat<- exposure_dat %>% add_column(samplesize.exposure=115078,
                                             units.exposure="SD")  
  exposure<-exposures[i]
  

#####################################################################################################################################
#2. Extract SNPs from outcome Data
######################################################################################################################################
#1. MDD outcome data: PGC/ 23andme outcome minus UKBB  
#2. rMDD outcome data: UKBB recurrent MDD sample 

    outcome_dat <-  read_outcome_data(
    filename= "./MDDsumstats.txt", 
    snps = exposure_dat$SNP, 
    snp_col = "SNP", 
    beta_col = "OR", 
    sep = "\t", 
    se_col = "SE", 
    pval_col = "P", 
    eaf_col =  "EAF", 
    effect_allele_col = "EA", 
    other_allele_col = "NEA"
  )
  
  outcome_dat$beta.outcome<- log(outcome_dat$beta.outcome)
  outcome_dat<- outcome_dat %>% add_column(samplesize.outcome=480359,
                                           ncase.outcome= 135458,units.outcome="log odds",
                                           ncontrol.outcome=344901, prevalence.outcome=0.1)  ### for PGC MDD 807k
  

  outcome<- "MDD807_minusUKBB"
  
######################################################################################################################################
#3. Harmonise data
######################################################################################################################################
  
  dat <- harmonise_data(exposure_dat , 
                        outcome_dat, 
                        action = 2)
  
  dat<- subset(dat, dat$mr_keep==T)
  dat<-clump_data(dat, clump_r2 = 0.001)
  
######################################################################################################################################
#4. Calculate F statistics
######################################################################################################################################
  
  dat$betaXG<-dat$beta.exposure
  dat$seBetaXG <-dat$se.exposure
  
  BetaXG=dat$betaXG
  seBetaXG = dat$seBetaXG
  BXG = abs(BetaXG)
  
  dat$F=dat$beta.exposure^2/dat$se.exposure^2

  
######################################################################################################################################
#5. Run MR
######################################################################################################################################
  
  mr_results <- mr(dat) 
  mr_results <- generate_odds_ratios(mr_results)

  mr_steiger<-steiger_filtering(dat)
  mr_raps<- as.data.frame(mr_raps(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome))
  
  mr_pleio<- mr_pleiotropy_test(dat)
  mr_hetero<- mr_heterogeneity(dat)

  mr_single <- mr_singlesnp(dat)
  mr_loo<- mr_leaveoneout(dat)

  
######################################################################################################################################
#6. Save results
######################################################################################################################################

  path= paste0("./Results/", exposure, "&", outcome, ".xlsx")
   
  
  write.xlsx(dat, file = path, sheetName = "dat")
  write.xlsx(mr_results, file = path, sheetName = "mr_results", append=T)
  write.xlsx(F, file = path, sheetName = "F", append=T)
  write.xlsx(mr_pleio, file = path, sheetName = "mr_pleio", append=T)
  write.xlsx(mr_hetero, file = path, sheetName = "mr_hetero", append=T)
  write.xlsx(mr_steiger, file = path, sheetName = "mr_steiger", append=T)
  write.xlsx(mr_loo, file = path, sheetName = "mr_loo", append=T)
  write.xlsx(mr_single, file = path, sheetName = "mr_single", append=T)
  write.xlsx(mr_raps, file = path, sheetName = "mr_raps", append=T)

  ####This html file contains plots for above analyses:
  path =  paste0("./Results/", exposure, "&", outcome,".html")
  mr_report(dat, output_path=path) 
    
  rm(list=ls(all=TRUE)) #empties your R environment
  }

