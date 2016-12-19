#1415 Col4a5xDO qtl 
#Yuka Takemon
#11/01/16

#This script will calculate the sex effect of Alb qtl map, by using sex as an interactive, while keeping creatinine as an additive covariate.

#load packages
library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load files
#genoprobs
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
#GM_grid
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
#kinship
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
#phenotype
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean up phenotype data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA

pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Alb10WK_log <- log(pheno$Alb10WK)
pheno$Alb15WK_log <- log(pheno$Alb15WK)

pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno$Creat10WK_log <- log(pheno$Creat10WK)
pheno$Creat15WK_log <- log(pheno$Creat15WK)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#create covar with sex and creatinine
#All 6wk
temp <- sex.covar
temp$creat6wk <- pheno$Creat6WK_log
covar.sex.creat6wk <- temp
#All 10wk
temp <- sex.covar
temp$creat10wk <- pheno$Creat10WK_log
covar.sex.creat10wk <- temp
#All 15wk
temp <- sex.covar
temp$creat15wk <- pheno$Creat15WK_log
covar.sex.creat15wk <- temp

#create qtl object
qtl.log.Alb6WK.192.sex.intcovar <- scanone(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat6wk, intcovar = sex.covar, snps = GM_snps)
qtl.log.Alb10WK.192.sex.intcovar <- scanone(pheno = pheno, pheno.col = "Alb10WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat10wk, intcovar = sex.covar, snps = GM_snps)
qtl.log.Alb15WK.192.sex.intcovar <- scanone(pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat15wk, intcovar = sex.covar, snps = GM_snps)

#save qtl 
save(qtl.log.Alb6WK.192.sex.intcovar, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.sex.intcovar.Rdata")
save(qtl.log.Alb10WK.192.sex.intcovar, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.sex.intcovar.Rdata")
save(qtl.log.Alb15WK.192.sex.intcovar, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.sex.intcovar.Rdata")

#plot each qlt
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.Alb6wk_log.sex.intcovar.png", width = 2000, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192.sex.intcovar, main = "Col4a5xDO logAlb6wk QTL map w/ creatinine as addcovar & sex as intcovar")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.Alb10wk_log.sex.intcovar.png", width = 2000, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192.sex.intcovar, main = "Col4a5xDO logAlb10wk QTL map w/ creatinine as addcovar & sex as intcovar")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.Alb15wk_log.sex.intcovar.png", width = 2000, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192.sex.intcovar, main = "Col4a5xDO logAlb15wk QTL map w/ creatinine as addcovar & sex as intcovar")
dev.off()

#Permutation for QTL with interactive covaritate is no working.
#For manyal permutaGo to: QTL_sex_intcovar_Alb_perm.R 

##################################
#1000.perms.txt was created with 1000 permutaitons with an interactive covariate. 

#load compiled perms data
perms.1000.qtl.Alb6.sexint <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.sex.int.perm/1000.perms.txt")
perms.1000.qtl.Alb10.sexint <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.sex.int.perm/1000.perms.txt")
perms.1000.qtl.Alb15.sexint <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.sex.int.perm/1000.perms.txt")

#perms must be in matrix form for get.sig.thr()
perms.1000.qtl.Alb6.sexint <- as.matrix(perms.1000.qtl.Alb6.sexint)
perms.1000.qtl.Alb10.sexint <- as.matrix(perms.1000.qtl.Alb10.sexint)
perms.1000.qtl.Alb15.sexint <- as.matrix(perms.1000.qtl.Alb15.sexint)

#get threshold
sig.thr_6 <- get.sig.thr(perms.1000.qtl.Alb6.sexint, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )
sig.thr_10 <- get.sig.thr(perms.1000.qtl.Alb10.sexint, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )
sig.thr_15 <- get.sig.thr(perms.1000.qtl.Alb15.sexint, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )

#plot again with threshold
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qlt.Alb6.sex.intcovar.perm1000.png", width = 2000, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192.sex.intcovar, sig.thr = sig.thr_6, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO logAlb6wk QTL map w/ creatinine as addcovar & sex as intcovar w/ 1000 perms")
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qlt.Alb10.sex.intcovar.perm1000.png", width = 2000, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192.sex.intcovar, sig.thr = sig.thr_10, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO logAlb10wk QTL map w/ creatinine as addcovar & sex as intcovar w/ 1000 perms")
dev.off()


png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qlt.Alb15.sex.intcovar.perm1000.png", width = 2000, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192.sex.intcovar, sig.thr = sig.thr_15, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO logAlb15wk QTL map w/ creatinine as addcovar & sex as intcovar w/ 1000 perms")
dev.off()

#having trouble plotting both qlt w/ and w/o sex interaction into one plot.... Asking dan 
#DOQTL does not handle plotting both maps

#try plotting Hi = Hf - Ha to show difference/interaction
#to perserve DOQTL class format, create a fake QTL object and replace just its LOD scores for mapping purposes. 
#not sure how to do this for its coefficient value yet... 

#The Hi is calculated by Hi = Hf - Ha, but only for the LOD scores
#All other values of Hi is "fake". Hi object is a shell which retains DOQTL formatting to hold unique LOD scores for compatibility


#load Hf and Ha qtl objects
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.Rdata")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.sex.intcovar.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.sex.intcovar.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.sex.intcovar.Rdata")

Hf_6 <- qtl.log.Alb6WK.192.sex.intcovar
Hf_10 <- qtl.log.Alb10WK.192.sex.intcovar
Hf_15 <- qtl.log.Alb15WK.192.sex.intcovar

Ha_6 <- qtl.log.ACR6WK.192
Ha_10 <- qtl.log.ACR10WK.192
Ha_15 <-  qtl.log.ACR15WK.192

#Assign shell for Hi 
Hi_6 <- Hf_6
Hi_10 <- Hf_10
Hi_15 <- Hf_15

#Calculate Hi = Hf - Ha
Hi_6$lod$A$lod <- Hi_6$lod$A$lod - Ha_6$lod$A$lod
Hi_10$lod$A$lod <- Hi_10$lod$A$lod - Ha_10$lod$A$lod
Hi_15$lod$A$lod <- Hi_15$lod$A$lod - Ha_15$lod$A$lod

#save Hi objects
save(Hi_6, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/Hi_qtl.log.Alb6WK.192.Rdata")
save(Hi_10, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/Hi_qtl.log.Alb10WK.192.Rdata")
save(Hi_15, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/Hi_qtl.log.Alb15WK.192.Rdata")

##run permutation to test for significance in: Hi_QTL_Alb_perm.R

#############################
#1000.perms.txt was created with 1000 permutaitons of Hi
Alb6.perm1000.Hi <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.Hi.perm/1000.perms.txt", header = FALSE)
Alb10.perm1000.Hi <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.Hi.perm/1000.perms.txt", header = FALSE)
Alb15.perm1000.Hi <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.Hi.perm/1000.perms.txt", header = FALSE)

#save 1000 perms
save(Alb6.perm1000.Hi, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.192.Hi.Rdata")
save(Alb10.perm1000.Hi, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.192.Hi.Rdata")
save(Alb15.perm1000.Hi, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15WK.192.Hi.Rdata")

#get sig thresholds
sig.thr_6 <- get.sig.thr(Alb6.perm1000.Hi, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )
sig.thr_10 <- get.sig.thr(Alb10.perm1000.Hi, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )
sig.thr_15 <- get.sig.thr(Alb15.perm1000.Hi, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )

#plot with significant thresholds
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.Alb6WK.Hi.perm1000.png", width = 2000, height = 1000, res = 100)
plot(Hi_6, sig.thr = sig.thr_6, sig.col = c("red", "orange", "chartreuse"), main = "Hi: Col4a5xDO logAlb6WK Interaction only QTL map w/ 1000 perms")
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.Alb10WK.Hi.perm1000.png", width = 2000, height = 1000, res = 100)
plot(Hi_10, sig.thr = sig.thr_10, sig.col = c("red", "orange", "chartreuse"), main = "Hi: Col4a5xDO logAlb10WK Interaction only QTL map w/ 1000 perms")
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.Alb15WK.Hi.perm1000.png", width = 2000, height = 1000, res = 100)
plot(Hi_15, sig.thr = sig.thr_15, sig.col = c("red", "orange", "chartreuse"), main = "Hi: Col4a5xDO logAlb15WK Interaction only QTL map w/ 1000 perms")
dev.off()












