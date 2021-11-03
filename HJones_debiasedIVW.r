# MRBase package
.libPaths("C:/R/Library")
# To update the package just run:
library(usethis)
library(devtools)
install_github("MRCIEU/TwoSampleMR")
install_github("tye27/mr.divw")
library(TwoSampleMR)
library(ggplot2)
library(mr.divw)

##############################################################################################################################################################################################
### PERFORMING DEBIASED INVERSE-VARIANCE WEIGHTED ESTIMATOR IN 2-SAMPLE MR
### REFERENCE: https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-4/Debiased-inverse-variance-weighted-estimator-in-two-sample-summary-data/10.1214/20-AOS2027.full
### GITHUB SCRIPT AND WIKI: https://github.com/tye27/mr.divw
##############################################################################################################################################################################################

#####################################################################
### 5e-5 threshold, DHA = 44 SNPS
#####################################################################
# This is a dataframe I've set up to collect results as I was performing this for lots of fatty acids at the time
rm(list=ls(all=TRUE))
results<-data.frame()

#####################################################################
#####################################################################
### DHA
#####################################################################
#####################################################################
# working directory
setwd("C:/Users/Documents")

# exposure data - Kettunen IVs, CHARGE effects, space delimited file with following headings:
# Phenotype SNP effect_allele other_allele eaf beta se pval
exp_dat <- read_exposure_data("DHA_IVs_DHA_effects_exposure_info.txt")
head(exp_dat)

# outcome data - space delimited file with following headings (note, not all of these are needed!):
# Phenotype SNP CHR BP effect_allele other_allele eaf beta se pval
out_dat <- read_outcome_data("SCZ_DHA_outcome_info_MRBaseclump.txt")
head(out_dat)

# list of SNPs with corresponding P values (taken from Kettunen statistics)
# space delimited file - 2 columns (SNP and p-value)
KetPs <- read.table("Kettunen_DHA_IVs_Pvalues.txt")
colnames(KetPs) <- c("SNP", "KetP")

#####################################################################
# dIVW MR method SCZ ~ DHA
#######################################################################
# harmonising data using TwoSampleMR package
dat <- harmonise_data( 
	exposure_dat = exp_dat,
	outcome_dat = out_dat
)

# merge in Kettunen p-values
dat1 <- merge(dat, KetPs, by="SNP")

# mr.divw
mr.dvw <- mr.divw(dat1$beta.exposure, dat1$beta.outcome, dat1$se.exposure, dat1$se.outcome, pval.selection=KetP)
z <- mr.dvw$beta.hat / mr.dvw$beta.se
p <- 2*pnorm((abs(z)), lower.tail=F)
mr.dvw
p

# Pushing results to dataframe (again, this was more useful when I was looking at multiple exposures)
results <-rbind(results, c("DHA",mr.dvw$n.IV, mr.dvw$beta.hat, mr.dvw$beta.se, p, mr.dvw$condition, mr.dvw$tau.square))

####################
# naming results dataframe columns
colnam <-c("FA","nsnp","b","se","pval","condition","tau.square")
colnames(results) <- colnam
# w3 results
results
