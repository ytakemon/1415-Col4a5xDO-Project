#This is to compare 3 differnt genoprobs files to see if our data is accurate
#Genoprobs 1: 157GS + 24RS + 11Gillian = 192 samples (best.genoprobs.192)
#Genoprobs 2: Only 157 GS genoprobs (genoprobs.GS.157)
#Genoprobs 3: 182GS (genoprobs.GS.182)

#source and laod DOQTL package: this only needs to happen once
#source("http://bioconductor.org/biocLite.R")
#install.packages("devtools")
#install_github("dmgatti/DOQTL")
#install_github("yihui/knitr")
#install_github("cran/abind")

#Load packages 
library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load 3 genoprobs
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GS.genoprobs.157.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GS.genoprobs.182.Rdata")

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno$C2[123] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2) 
pheno$X2xC1_log <- log(pheno$X2xC1)
pheno$C1_log <- log(pheno$C1)
pheno$PL_log <- log(pheno$PL)
options(na.action = 'na.pass') #leave in NAs

#create covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#submit following sections separately to compare 3 datasets
#QTL mapping
#genoprob 192
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
qtl.GFR.log.C2.192 <- scanone( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.GFR.log.C2.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C2.192.Rdata" )

#genoprob 157
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(GS.genoprobs.157, snps = GM_snps, bychr = TRUE)
qtl.GFR.log.C2.157 <- scanone( pheno = pheno, pheno.col = "C2_log", probs = GS.genoprobs.157, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.GFR.log.C2.157, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C2.157.Rdata" )

#genoprob 182
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(GS.genoprobs.182, snps = GM_snps, bychr = TRUE)
qtl.GFR.log.C2.182 <- scanone( pheno = pheno, pheno.col = "C2_log", probs = GS.genoprobs.182, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.GFR.log.C2.182, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C2.182.Rdata" )


########################################################################################
#plots qtl on Rstudio
tiff("./genorpob192.logC2GFR.qtl.tiff")
plot(qtl.GFR.log.C2.192, main ="Col4a5xDO genoprob.192 logC2.GFR QTL")
dev.off()
tiff("./genorpob157.logC2GFR.qtl.tiff")
plot(qtl.GFR.log.C2.157, main ="Col4a5xDO genoprob.157 logC2.GFR QTL")
dev.off()
tiff("./genorpob182.logC2GFR.qtl.tiff")
plot(qtl.GFR.log.C2.182, main ="Col4a5xDO genoprob.182 logC2.GFR QTL")
dev.off()
########################################################################################

#QTL perms 100
#genoprob 192
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
perms.100.qtl.GFR.log.C2.192 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 100)
save(perms.100.qtl.GFR.log.C2.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.C2.192.Rdata" ) 

#genoprob 157
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(GS.genoprobs.157, snps = GM_snps, bychr = TRUE)
perms.100.qtl.GFR.log.C2.157 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = GS.genoprobs.157, addcovar = sex.covar, snps = GM_snps, nperm = 100)
save(perms.100.qtl.GFR.log.C2.157, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.C2.157.Rdata" ) 

#genoprob 182
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
perms.100.qtl.GFR.log.C2.182 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = GS.genoprobs.182, addcovar = sex.covar, snps = GM_snps, nperm = 100)
save(perms.100.qtl.GFR.log.C2.182, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.C2.182.Rdata" ) 


#QTL perms 1000
#genoprob 192
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
perms.1000.qtl.GFR.log.C2.192 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, addcovar = sex_covar, snps = GM_snps, nperm = 100)
save(perms.1000.qtl.GFR.log.C2.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perms1000/perms.1000.qtl.GFR.log.C2.192.Rdata" ) 

#genoprob 157
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(GS.genoprobs.157, snps = GM_snps, bychr = TRUE)
perms.1000.qtl.GFR.log.C2.157 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = GS.genoprobs.157, addcovar = sex_covar, snps = GM_snps, nperm = 100)
save(perms.1000.qtl.GFR.log.C2.157, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perms1000/perms.1000.qtl.GFR.log.C2.192.Rdata" ) 

#genoprob 182
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
perms.1000.qtl.GFR.log.C2.182 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = GS.genoprobs.182, addcovar = sex_covar, snps = GM_snps, nperm = 100)
save(perms.1000.qtl.GFR.log.C2.182, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perms1000/perms.1000.qtl.GFR.log.C2.182.Rdata" ) 

#with perms create thresholds
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perms.100.qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perms.100.qtl.GFR.log.C2.157.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perms.100.qtl.GFR.log.C2.182.Rdata")

thr.100.qtl.GFR.logC2.192 <- get.sig.thr( perms.100.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.100.qtl.GFR.logC2.157 <- get.sig.thr( perms.100.qtl.GFR.log.C2.157[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.100.qtl.GFR.logC2.182 <- get.sig.thr( perms.100.qtl.GFR.log.C2.182[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

save(thr.100.qtl.GFR.logC2.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/thr.100.qtl.GFR.logC2.192.Rdata")
save(thr.100.qtl.GFR.logC2.157, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/thr.100.qtl.GFR.logC2.157.Rdata")
save(thr.100.qtl.GFR.logC2.182, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/thr.100.qtl.GFR.logC2.182.Rdata")

##########################################################################################################################
#on Rstudio
tiff("./genoprob192.logC2GFR.qtl.tiff")
plot(qtl.GFR.log.C2.192, sig.thr = thr.100.qtl.GFR.logC2.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 192 genoprob log.C2.GFR QTL perms_100")
dev.off()

tiff("./genoprob182.logC2GFR.qtl.tiff")
plot(qtl.GFR.log.C2.182, sig.thr = thr.100.qtl.GFR.logC2.182, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 182 genoprob log.C2.GFR QTL perms_100")
dev.off()

tiff("./genoprob157.logC2GFR.qtl.tiff")
plot(qtl.GFR.log.C2.157, sig.thr = thr.100.qtl.GFR.logC2.157, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 157 genoprob log.C2.GFR QTL perms_100")
dev.off()


#conclusion: THANK GOD WE DID THIS! Chr 7 was definetly a false peak, and peak on Chr 19 was masked by other bad quality data.

