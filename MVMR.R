######################################################################################################################
# Multivariable MR Omega 3 Fatty Acids
######################################################################################################################

rm(list=ls(all=TRUE)) #empties your R environment
ao<-available_outcomes()
setwd("/")

source("/R/Source code/InstallPackages.R") # installs required R packages- also given in Omega 3 MR.R
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
library(TwoSampleMR)

######################################################################################################################################
# 1. Call Exposure Data
######################################################################################################################################
######################################################################################################################################

######################################################################################################################################
# MVMR MODEL 1 (Omega 3 & 6)
######################################################################################################################################

exposure_dat<-mv_extract_exposures(id_exposure = c("met-d-Omega_3", "met-d-Omega_6")) 

######################################################################################################################################
# MVMR MODEL 2 (Omega 3, HDL, LDL, Trigs)
######################################################################################################################################

exposure_dat<-mv_extract_exposures(id_exposure = c("met-d-Omega_3", "ieu-b-111", "ieu-b-109", "ieu-b-110")) 

######################################################################################################################################
# MVMR MODEL 3 (EPA & DHA)
######################################################################################################################################

## Model 3 (ie EPA/ DHA) uses UKBB SNPs with CHARGE EPA and DHA SNP-exposure estimates (as EPA is not measured in UKBB)

## a. Using UK Biobank SNPs:

### i. Extract UK Biobank SNPs

exposure_dat<- extract_instruments(outcomes= "met-d-Omega_3", clump = T)
snps <- exposure_dat$SNP

### ii. extract SNP-exposure data from EPA and DHA CHARGE files

EPA<- TwoSampleMR::read_outcome_data(filename = "./CHARGE_N3_EPA_fixed.txt",snps = snps,
                                     sep = "\t", beta_col = "beta.sd", se_col = "se.sd", pval_col = "pval.sd")
EPA<-convert_outcome_to_exposure(EPA)
EPA$exposure<- "EPA"

DHA<- TwoSampleMR::read_outcome_data(filename = "./CHARGE_N3_DHA_fixed.txt", snps=snps,
                                     sep = "\t", beta_col = "beta.sd", se_col = "se.sd", pval_col = "pval.sd")
DHA<- convert_outcome_to_exposure(DHA)
DHA$exposure<- "DHA"

exposure_dat<-rbind(DHA, EPA)

##b. Using Kettunen SNPs (as these improve conditional F statistics)
# i. Extract Kettunen Omega 3 SNPs (as above but using kettunen GWAS)
exposure_dat<-extract_instruments("met-c-855", clump=F)

## ii. as above for extracting from EPA/ DHA 

######################################################################################################################################
#2. Call Outcome Data
######################################################################################################################################
######################################################################################################################################
# As per Omega 3 MR.R
#source("./Extract SNPs MDD.R")

source("./Extract SNPs rMDD.R")

######################################################################################################################################
#3. Format and Harmonise data
######################################################################################################################################
######################################################################################################################################

dat<-mv_harmonise_data(exposure_dat, outcome_dat,2)

######################################################################################################################################
#4. Results
######################################################################################################################################
######################################################################################################################################

res <- mv_multiple(dat, plots = T)
res
#MVMR<- as.data.frame(res)
#MVMR

write.xlsx(MVMR, "./MVMR.xlsx", sheetName = "model123", append=T)

######################################################################################################################################
#5. Conditional F stats
######################################################################################################################################
######################################################################################################################################
# MODEL 1
######################################################################################################################################

n3 <- subset(exposure_dat, exposure_dat$id.exposure=="met-d-Omega_3")
n3<- n3[order(n3$SNP),]
n3<-subset(n3, n3$SNP %in% outcome_dat$SNP)

n6 <- subset(exposure_dat, exposure_dat$id.exposure=="met-d-Omega_6")
n6<- n6[order(n6$SNP),]
n6<-subset(n6, n6$SNP %in% outcome_dat$SNP)

n3<- n3[order(n3$SNP),]
n6<- n6[order(n6$SNP),]

mvmrdat<-cbind(n3$SNP, outcome_dat$beta.outcome, outcome_dat$se.outcome,
               n3$beta.exposure, n3$se.exposure,
               n6$beta.exposure, n6$se.exposure)


F.dat<- format_mvmr(BXGs = mvmrdat[,c(4,6)],
                    BYG = mvmrdat[,2],
                    seBXGs = mvmrdat[,c(5,7)],
                    seBYG = mvmrdat[,3],
                    RSID = mvmrdat[,1])

bX1<- c(dat$exposure_beta[,1])
bX2<- c(dat$exposure_beta[,2])
bY<- c(dat$outcome_beta)
bXse1<- c(dat$exposure_se[,1])
bXse2<- c(dat$exposure_se[,2])
bYse<- c(dat$outcome_se)

# CORRELATIONS
# n3 x n6 = 0.77

#     n3 n6
# n3 1   0.77
# n6 0.77 1

# Covariance between exposures
corrmat <- matrix(, nrow = 2, ncol = 2)

colnames(corrmat) <- c("n3","n6")
rownames(corrmat) <- c("n3","n6")

corrmat[c("n3", "n6"),c("n6","n3")] <- 0.77

diag(corrmat) <- 1

Xse <- cbind(bXse1,bXse2)

colnames(Xse) <- c("n3", "n6")

corrmat
head(Xse)

df<- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)
i<- format_mvmr(df[,c(1,3)],
                df[,5],
                df[,c(2,4)],
                df[,6])


# note that phenocov changes correlation matrix into covariance matrix
cov <-phenocov_mvmr(corrmat, Xse) # <- may need to change to covariance matrix instead...

x<- strength_mvmr(i, cov) # the strength mvmr command will give you the F statistic

z<- MVMR::ivw_mvmr(r_input = F.dat, gencov = cov)

######################################################################################################################################
# Model 2 (Omega 3, HDL, LDL, Trigs)
######################################################################################################################################

exposure_dat<- subset(exposure_dat, exposure_dat$SNP %in% outcome_dat$SNP)
  
exposure_dat<-exposure_dat[order(exposure_dat$SNP),]
outcome_dat<-outcome_dat[order(outcome_dat$SNP),]
     
#rmSNP<- c("rs10822163", "rs11722924", "rs5112","rs551243", "rs5755799", "rs62112763", "rs6752845", "rs6916318", "rs9368503") ## removed as palindromic
#SNPs <- setdiff(SNPs,rmSNP) 

n3 <- subset(exposure_dat, exposure_dat$id.exposure=="met-d-Omega_3")
n3<- n3[order(n3$SNP),]
n3<-subset(n3, n3$SNP %in% outcome_dat$SNP)

HDL <- subset(exposure_dat, exposure_dat$id.exposure=="ieu-b-109")
HDL<- HDL[order(HDL$SNP),]
HDL<-subset(HDL, HDL$SNP %in% outcome_dat$SNP)

LDL <- subset(exposure_dat, exposure_dat$id.exposure=="ieu-b-110")
LDL<- LDL[order(LDL$SNP),]
LDL<-subset(LDL, LDL$SNP %in% outcome_dat$SNP)

Trig <- subset(exposure_dat, exposure_dat$id.exposure=="ieu-b-111")
Trig<- Trig[order(Trig$SNP),]
Trig<-subset(Trig, Trig$SNP %in% outcome_dat$SNP)


mvmrdat<-cbind(n3$SNP, outcome_dat$beta.outcome, outcome_dat$se.outcome,
               n3$beta.exposure, n3$se.exposure,
               Trig$beta.exposure,Trig$se.exposure,
               HDL$beta.exposure, HDL$se.exposure,
               LDL$beta.exposure, LDL$se.exposure)

F.dat<- format_mvmr(BXGs = mvmrdat[,c(4,6,8,10)],
                    BYG = mvmrdat[,2],
                    seBXGs = mvmrdat[,c(5,7,9,11)],
                    seBYG = mvmrdat[,3],
                    RSID = mvmrdat[,1])

#		HDL		LDL		TG		Omega3
# HDL	1		0.1		-0.437	-0.0667
# LDL	0.1		1		0.222	0.49
# TG	-0.437	0.222	1		0.58
# Omega3	-0.667	0.49	0.58	1


bX1<- c(dat$exposure_beta[,1])
bX2<- c(dat$exposure_beta[,2])
bX3<- c(dat$exposure_beta[,3])
bX4<- c(dat$exposure_beta[,4])
bY<- c(dat$outcome_beta)

bXse1<- c(dat$exposure_se[,1])
bXse2<- c(dat$exposure_se[,2])
bXse3<- c(dat$exposure_se[,3])
bXse4<- c(dat$exposure_se[,4])
bYse<- c(dat$outcome_se)


df<- data.frame(bX1, bXse1, bX2, bXse2, bX3, bXse3, bX4, bXse4, bY, bYse)
i<- format_mvmr(df[,c(1,3,5,7)],df[,9],df[,c(2,4,6,8)],df[,10])

# Covariance between exposures
corrmat <- matrix(nrow = 4, ncol = 4)

colnames(corrmat) <- c("n3","Trig", "HDL", "LDL")
rownames(corrmat) <- c("n3","Trig", "HDL", "LDL")

corrmat[c("n3", "Trig"),c("n3","Trig")] <- 0.58 
corrmat[c("n3", "HDL"),c("n3","HDL")] <- -0.0668 
corrmat[c("n3", "LDL"),c("n3","LDL")] <- 0.49 
corrmat[c("HDL", "LDL"),c("LDL","HDL")] <- 0.1
corrmat[c("HDL", "Trig"),c("HDL","Trig")] <- -0.437
corrmat[c("Trig", "LDL"),c("LDL","Trig")] <- 0.222
corrmat

diag(corrmat) <- 1

Xse <- cbind(bXse1,bXse2,bXse3,bXse4)

colnames(Xse) <- c("n3","Trig", "HDL", "LDL")
head(Xse)

# note that phenocov changes correlation matrix into covariance matrix
cov <-phenocov_mvmr(corrmat, Xse) # <- may need to change to covariance matrix instead...
cov
x<- strength_mvmr(i, cov) # the strength mvmr command will give you the F statistic

x<- strength_mvmr(i, cov) # the

#z<-mvmr(r_input = F.dat)

pleiotropy_mvmr(r_input = F.dat, gencov = 0)

######################################################################################################################################
# DHA/ EPA

#DHA <- subset(exposure_dat, exposure_dat$exposure=="DHAsd")
#DHA<- DHA[order(DHA$SNP),]
DHA<-subset(DHA, DHA$SNP %in% SNPs1)

#EPA <- subset(exposure_dat, exposure_dat$exposure=="EPAsd")
#EPA<- EPA[order(EPA$SNP),]
#EPA<-subset(EPA, EPA$SNP %in% SNPs)


mvmrdat<-cbind(exposure_dat$SNP, outcome_dat$beta.outcome, outcome_dat$se.outcome,
               EPA$beta.exposure, EPA$se.exposure,
               DHA$beta.exposure, DHA$se.exposure)

F.dat<- format_mvmr(BXGs = mvmrdat[,c(4,6)],
                    BYG = mvmrdat[,2],
                    seBXGs = mvmrdat[,c(5,7)],
                    seBYG = mvmrdat[,3],
                    RSID = mvmrdat[,1])
F.dat$SNP<-(1:82)
F.dat

# CORRELATIONS
# EPA x DHA = 0.460143966613421

#     EPA DHA
# EPA 1   0.46
# DHA 0.46 1

bX1<- c(dat$exposure_beta[,1])
bX2<- c(dat$exposure_beta[,2])
bY<- c(dat$outcome_beta)
bXse1<- c(dat$exposure_se[,1])
bXse2<- c(dat$exposure_se[,2])
bYse<- c(dat$outcome_se)

df<- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)
i<- format_mvmr(df[,c(1,3)],
                df[,5],
                df[,c(2,4)],
                df[,6])

corrmat <- matrix(nrow = 2, ncol = 2)

colnames(corrmat) <- c("EPA","DHA")
rownames(corrmat) <- c("EPA","DHA")

corrmat[c("EPA", "DHA"),c("DHA","EPA")] <- 0.46 
diag(corrmat) <- 1

Xse <- cbind(bXse1,bXse2)

colnames(Xse) <- c("EPA","DHA")
head(Xse)

# note that phenocov changes correlation matrix into covariance matrix
cov <-phenocov_mvmr(corrmat, Xse) # <- may need to change to covariance matrix instead...
cov
x<- strength_mvmr(i, cov) # the strength mvmr command will give you the F statistic

z<-mvmr(i)

pleiotropy_mvmr(r_input =i, gencov = cov)


######################################################################################################################################







