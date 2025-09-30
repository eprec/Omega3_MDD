#######################################################################################################################
## OMEGA 3 MDD MR SCRIPT
#######################################################################################################################

#######################################################################################################################
## Set up
#######################################################################################################################

rm(list=ls(all=TRUE)) #empties your R environment
setwd("/")
source("R/Source code/InstallPackages.R")

install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
library(TwoSampleMR)

#######################################################################################################################
##1.  Create exposure data
#######################################################################################################################
# call UKBB SNPs
exposure_dat<- extract_instruments(outcomes= "met-d-Omega_3", clump = F)

#exposure_dat<- extract_instruments(outcomes= "met-c-855", clump = F)

snps <- exposure_dat$SNP

#extract these from CHARGE files
EPA<- TwoSampleMR::read_outcome_data(filename = "./CHARGE_N3_EPA_fixed.txt",snps = snps,
                                     sep = "\t", beta_col = "beta.sd", se_col = "se.sd", pval_col = "pval.sd")
EPA<-convert_outcome_to_exposure(EPA)
EPA$exposure<- "EPA using UKBB n3 SNPs"

exposure_dat<-EPA

#######################################################################################################################
##2.  Call outcome data
#######################################################################################################################
# As in Omega 3 MR.R

# label outcome and exposure for paths later

outcome<- outcome_dat$outcome[1]
exposure<- exposure_dat$exposure[1]

######################################################################################################################################
##3. Harmonise data
######################################################################################################################################

dat <- harmonise_data(exposure_dat , outcome_dat, action = 2)
dat<- subset(dat, dat$mr_keep==T)

dat<-clump_data(dat, clump_r2 = 0.001)

dat$r.exposure<- get_r_from_bsen(b = dat$beta.exposure, se = dat$se.exposure, n = 8866) 
dat$r.outcome<- get_r_from_lor(lor = dat$beta.outcome, af = dat$eaf.outcome, ncase = dat$ncase.outcome, ncontrol = dat$ncontrol.outcome, prevalence = 0.1, model = "logit")


######################################################################################################################################
## 4. F stats
######################################################################################################################################

dat$F=abs(dat$beta.exposure)^2/dat$se.exposure^2
mFc=mean(F)
sFc=sum(F)

######################################################################################################################################
## 5. Run MR
######################################################################################################################################

mr_results <- mr(dat) 
mr_results

mr_results <-subset(mr_results, mr_results$method != "Simple mode")
mr_results <- generate_odds_ratios(mr_results)

mr_results<- plyr::arrange(mr_results, mr_results$method)

mr_results$result <- paste(round(mr_results$or, digits = 2), 
                           "(", round(mr_results$or_lci95, digits = 2), 
                           "-", round(mr_results$or_uci95, digits = 2), ")") 

mr_results$pval <- signif(mr_results$pval, digits = 3)

mr_steiger<-steiger_filtering(dat)

mr_raps<- as.data.frame(mr_raps(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome))

mr_pleio<- mr_pleiotropy_test(dat)
mr_hetero<- mr_heterogeneity(dat)
mr_hetero<- plyr::arrange(mr_hetero, mr_hetero$method)
I2<- TwoSampleMR::Isq(y= dat$beta.exposure, s=dat$se.exposure)

mr_single <- mr_singlesnp(dat)
mr_loo<- mr_leaveoneout(dat)
mr_leaveoneout_plot(mr_single)

snpslost <- as.data.frame("nil")
snpslost <- as.data.frame(setdiff(exposure_dat$SNP, outcome_dat$SNP))



