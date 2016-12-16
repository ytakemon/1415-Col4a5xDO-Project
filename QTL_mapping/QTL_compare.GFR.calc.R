#comparing different GFR calculations to determine best calculation to use
#best one will be determined by making sure that they all show simialr data and have the most number of retained values
#eg. currently we have the most confidenence in the C2 calculation however due to its stringent QC, many samples are missing GFR values
library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load sample
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/best.genoprobs.192.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2) 
pheno$X2xC1_log <- log(pheno$X2xC1)
pheno$C1_log <- log(pheno$C1)
pheno$PL_log <- log(pheno$PL)
options(na.action = 'na.pass') #leave in NAs

#Check out distribution of GFR
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/hist.GFR.C2.compare.png", width = 1000, height = 1000, res = 100) 
layout(matrix(1:2,1,2))
hist(pheno$C2)
hist(pheno$C2_log)
dev.off()


#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"


#submit following sections separately to compare the 4 different GFR calculations
#QTL mapping
#logC2 GFR is already done
qtl.GFR.log.C2.192 <- scanone( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#logX2xC1 GFR
qtl.GFR.log.2xC1.192 <- scanone( pheno = pheno, pheno.col = "X2xC1_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#logC1 GFR
qtl.GFR.log.C1.192 <- scanone( pheno = pheno, pheno.col = "C1_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#logPL GFR
qtl.GFR.log.PL.192 <- scanone( pheno = pheno, pheno.col = "PL_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#save
save(qtl.GFR.log.C2.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C2.192.Rdata" )
save(qtl.GFR.log.2xC1.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.2xC1.192.Rdata" )
save(qtl.GFR.log.C1.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C1.192.Rdata" )
save(qtl.GFR.log.PL.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.PL.192.Rdata" )
#QTL mapping of female and male subset
#Female
qtl.GFR.log.C2.F <- scanone( pheno = pheno.F, pheno.col = "C2_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
qtl.GFR.log.2xC1.F <- scanone( pheno = pheno.F, pheno.col = "X2xC1_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
#Male
qtl.GFR.log.C2.M <- scanone( pheno = pheno.M, pheno.col = "C2_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
qtl.GFR.log.2xC1.M <- scanone( pheno = pheno.M, pheno.col = "X2xC1_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
#save
save(qtl.GFR.log.C2.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C2.F.Rdata")
save(qtl.GFR.log.2xC1.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.2xC1.F.Rdata")
save(qtl.GFR.log.C2.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.C2.M.Rdata")
save(qtl.GFR.log.2xC1.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.GFR.log.2xC1.M.Rdata")
#QTL mapping for C1_log of samples only in C2 pheno
qtl.GFR.log.C1.subset <- scanone( pheno = X2xC1, pheno.col = "X2xC1_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
qtl.GFR.log.C1.subset.F <- scanone( pheno = X2xC1, pheno.col = "X2xC1_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
qtl.GFR.log.C1.subset.M <- scanone( pheno = X2xC1, pheno.col = "X2xC1_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
#save
save(qtl.GFR.log.C1.subset, file =  "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C1.subset.Rdata")
save(qtl.GFR.log.C1.subset.F, file =  "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C1.subset.F.Rdata" )
save(qtl.GFR.log.C1.subset.M, file =  "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C1.subset.M.Rdata" )

#submit following sections separately to compare the 4 different GFR calculations
#run permutation for 1000
#C2 GFR
perms.1000.qtl.GFR.log.C2.192 <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
#2xC1 GFR
perms.1000.qtl.GFR.log.2xC1.192 <- scanone.perm( pheno = pheno, pheno.col = "X2xC1_log", probs = best.genoprobs.192,addcovar = sex.covar, snps = GM_snps, nperm = 1000)
#C1 GFR
perms.1000.qtl.GFR.log.C1.192 <- scanone.perm( pheno = pheno, pheno.col = "C1_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
#PL GFR
perms.1000.qtl.GFR.log.PL.192 <- scanone.perm( pheno = pheno, pheno.col = "PL_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
#save
save(perms.1000.qtl.GFR.log.C2.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C2.192.Rdata" ) 
save(perms.1000.qtl.GFR.log.2xC1.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.2xC1.192.Rdata" ) 
save(perms.1000.qtl.GFR.log.C1.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C1.192.Rdata" ) 
save(perms.1000.qtl.GFR.log.PL.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.PL.192.Rdata" ) 

##########
#plot
#load qtl data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.2xC1.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C1.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.PL.192.Rdata")
#load perms 100
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.2xC1.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.C1.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.GFR.log.PL.192.Rdata")
#load perms 1000
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.2xC1.192.Rdata")
#thresholding 100
thr.100.qtl.GFR.log.C2.192 <- get.sig.thr( perms.100.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.100.qtl.GFR.log.2xC1.192 <- get.sig.thr( perms.100.qtl.GFR.log.2xC1.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.100.qtl.GFR.log.C1.192 <- get.sig.thr( perms.100.qtl.GFR.log.C1.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.100.qtl.GFR.log.PL.192 <- get.sig.thr( perms.100.qtl.GFR.log.PL.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#thresholding 1000
thr.1000.qtl.GFR.log.C2.192 <- get.sig.thr( perms.1000.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.GFR.log.2xC1.192 <- get.sig.thr( perms.1000.qtl.GFR.log.2xC1.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#Plot 100
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perms100.GFR.log.C2.192.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C2.192, sig.thr = thr.100.qtl.GFR.log.C2.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO perms100 192 genoprob log.C2.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perms100.GFR.log.2xC1.192.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.2xC1.192, sig.thr = thr.100.qtl.GFR.log.2xC1.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO perms100 192 genoprob log.2xC1.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perms100.GFR.log.C1.192.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C1.192, sig.thr = thr.100.qtl.GFR.log.C1.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO perms100 192 genoprob log.C1.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perms100.GFR.log.PL.192.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.PL.192, sig.thr = thr.100.qtl.GFR.log.PL.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO perms100 192 genoprob log.PL.GFR QTL")
dev.off()
#Plot 1000
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms1000.GFR.log.C2.192.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C2.192, sig.thr = thr.1000.qtl.GFR.log.C2.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO perms1000 192 genoprob log.C2.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms1000.GFR.log.2xC1.192.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.2xC1.192, sig.thr = thr.1000.qtl.GFR.log.2xC1.192, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO perms1000 192 genoprob log.2xC1.GFR QTL")
dev.off()
#for female and male subset
#female
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.C2.F.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C2.F, main = "Col4a5xDO 64 Female log.C2.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.2xC1.F.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.2xC1.F, main = "Col4a5xDO 60 Female log.2xC1.GFR QTL")
dev.off()
#male
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.C2.M.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C2.M, main = "Col4a5xDO 95 Male log.C2.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.2xC1.M.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.2xC1.M, main = "Col4a5xDO 93 Male log.2xC1.GFR QTL")
dev.off()
#for X2xC1 subseted to C2 samples
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.C1.subset.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C1.subset, main = "Col4a5xDO 124 samples log.2xC1.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.C1.subset.F.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C1.subset.F, main = "Col4a5xDO 65 female samples log.2xC1.GFR QTL")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.log.C1.subset.M.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C1.subset.M, main = "Col4a5xDO 60 male samples log.2xC1.GFR QTL")
dev.off()

####################
#scatter plot to compare GFR and ACR
#subset pheno for 124 samples
options(na.action = 'na.exclude')
pheno.124 <- pheno[rownames(C2),]
pheno.64.F <- pheno.124[pheno.124$Sex == "F",]
pheno.60.M <- pheno.124[pheno.124$Sex == "M",]
#GFR: C2
#GFR: C2 vs ACR 6wk
#both sex
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.124.ACR6wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.124$C2, y = pheno.124$ACR6WK, type = "p", main = "Col4a5xDO 124 C2 GFR cor 6wk ACR", xlab = "C2 GFR", ylab = "6wk ACR")
abline(fit <- lm(ACR6WK ~ C2, data = pheno.124), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))

legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$, digits =4)))
dev.off()
#females
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.64.F.ACR6wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.64.F$C2, y = pheno.64.F$ACR6WK, type = "p", main = "Col4a5xDO 64 Female C2 GFR cor 6wk ACR", xlab = "C2 GFR", ylab = "6wk ACR")
abline(fit <- lm(ACR6WK ~ C2, data = pheno.64.F), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#males
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.60.M.ACR6wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.60.M$C2, y = pheno.60.M$ACR6WK, type = "p", main = "Col4a5xDO 60 Male C2 GFR cor 6wk ACR", xlab = "C2 GFR", ylab = "6wk ACR")
abline(fit <- lm(ACR6WK ~ C2, data = pheno.60.M), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#GFR: C2 vs ACR 10wk
#both sex
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.124.ACR10wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.124$C2, y = pheno.124$ACR10WK, type = "p", main = "Col4a5xDO 124 C2 GFR cor 10wk ACR", xlab = "C2 GFR", ylab = "10wk ACR")
abline(fit <- lm(ACR10WK ~ C2, data = pheno.124), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#females
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.64.F.ACR10wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.64.F$C2, y = pheno.64.F$ACR10WK, type = "p", main = "Col4a5xDO 64 Female C2 GFR cor 10wk ACR", xlab = "C2 GFR", ylab = "10wk ACR")
abline(fit <- lm(ACR10WK ~ C2, data = pheno.64.F), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#males
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.60.M.ACR10wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.60.M$C2, y = pheno.60.M$ACR10WK, type = "p", main = "Col4a5xDO 60 Male C2 GFR cor 10wk ACR", xlab = "C2 GFR", ylab = "10wk ACR")
abline(fit <- lm(ACR10WK ~ C2, data = pheno.60.M), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#GFR: C2 vs ACR 15wk
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.124.ACR15wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.124$C2, y = pheno.124$ACR15WK, type = "p", main = "Col4a5xDO 124 C2 GFR cor 15wk ACR", xlab = "C2 GFR", ylab = "15wk ACR")
abline(fit <- lm(ACR15WK ~ C2, data = pheno.124), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#females
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.64.F.ACR15wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.64.F$C2, y = pheno.64.F$ACR15WK, type = "p", main = "Col4a5xDO 64 Female C2 GFR cor 15wk ACR", xlab = "C2 GFR", ylab = "15wk ACR")
abline(fit <- lm(ACR15WK ~ C2, data = pheno.64.F), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#males
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2.60.M.ACR15wk.png", width = 1000, height = 1000, res = 100)
plot(x = pheno.60.M$C2, y = pheno.60.M$ACR15WK, type = "p", main = "Col4a5xDO 60 Male C2 GFR cor 15wk ACR", xlab = "C2 GFR", ylab = "15wk ACR")
abline(fit <- lm(ACR15WK ~ C2, data = pheno.60.M), col = "red")
legend("topright", bty= "n", legend = paste("R^2 = ", format(summary(fit)$r.squared, digits =4)))
dev.off()
#GFR: 2xC1
#GFR: 2xC1 vs ACR 6wk
#GFR: 2xC1 vs ACR 10wk
#GFR: 2xC1 vs ACR 15wk


## Plotting QLT map as .pdf without X chromosome:
## Yuka Takemon
## 12/05/16

## Setup
library(DOQTL)
library(ggplot2)
library(reshape2)
setwd("/hpcdata/ytakemon/Col4a5xDO")
#	load files
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C2.192.Rdata")
thr.1000.qtl.GFR <- get.sig.thr( perms.1000.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

## GFR_C2 QTL plot
#	remove X chr from qtl object
qtl <- qtl.GFR.log.C2.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.GFR_C2.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.GFR, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.GFR_C2 QTL perms.1000.noX")
dev.off()













