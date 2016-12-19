#1415 Col4a5xDO qtl 
#Yuka Takemon
#11/01/16

#This script will calculate the sex effect of GFR C2 qtl map, by using sex as an interactive

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
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2) 
options(na.action = 'na.pass') #leave in NAs

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

qtl.GFR.log.C2.sex.intcovar.192 <- scanone(pheno = pheno, pheno.col = "C2_log", probs =  best.genoprobs.192, K = K_GS, addcovar = sex.covar, intcovar = sex.covar, snps = GM_snps)
save(qtl.GFR.log.C2.sex.intcovar.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.sex.intcovar.192.Rdata")

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.C2.sex.intcovar.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C2.sex.intcovar.192, main = "Col4a5xDO logC2GFR sex intcovar QTL map")
dev.off()

#Permutation for QTL with interactive covaritate is no working.
#For manyal permutaGo to: QTL_sex_intcovar_GFR_perm.R 

##################################
#1000.perms.txt was created with 1000 permutaitons with an interactive covariate. 

#load compiled perms data
perms.1000.qtl.GFR.C2.192.sexint <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.sex.int.perm/1000.perms.txt", header = FALSE)
#perms must be in matrix for get.sig.thr()
perms.1000.qtl.GFR.C2.192.sexint <- as.matrix(perms.1000.qtl.GFR.C2.192.sexint)
save(perms.1000.qtl.GFR.C2.192.sexint, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C2.192.sexint.Rdata")
#get threshold
sig.thr <- get.sig.thr(perms.1000.qtl.GFR.C2.192.sexint, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )

#plot again with threshold
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.C2.sex.intcovar.perm1000.png", width = 2000, height = 1000, res = 100)
plot(qtl.GFR.log.C2.sex.intcovar.192, sig.thr = sig.thr, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO logC2GFR sex intcovar QTL map w/ 1000 perms")
dev.off()

#having trouble plotting both qlt w/ and w/o sex interaction into one plot.... Asking dan 
#DOQTL does not handle plotting both maps

#try plotting Hi = Hf - Ha to show difference/interaction
#to perserve DOQTL class format, create a fake QTL object and replace just its LOD scores for mapping purposes. 
#not sure how to do this for its coefficient value yet... 

#The Hi is calculated by Hi = Hf - Ha, but only for the LOD scores
#All other values of Hi is "fake". Hi object is a shell which retains DOQTL formatting to hold unique LOD scores for compatibility


#load Hf and Ha qtl objects
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.sex.intcovar.192.Rdata") 
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata") 

Hf <- qtl.GFR.log.C2.sex.intcovar.192 
Ha <- qtl.GFR.log.C2.192
Hi <- Hf  #creates shell for Hi using Hf values

Hi$lod$A$lod <- Hi$lod$A$lod - Ha$lod$A$lod #replaces LOD score with LODf - LODa
save(Hi, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/Hi_qtl.GFR.log.C2.192.Rdata")


#run permutation to test for significance in: Hi_QTL_GFR_perm.R

#############################
#1000.perms.txt was created with 1000 permutaitons of Hi

#load compiled perms data
perms.1000.qtl.GFR.C2.192.Hi <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.Hi.perm/1000.perms.txt", header = FALSE)
#perms must be in matrix for get.sig.thr()
perms.1000.qtl.GFR.C2.192.Hi  <- as.matrix(perms.1000.qtl.GFR.C2.192.Hi )
save(perms.1000.qtl.GFR.C2.192.Hi , file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C2.192.Hi.Rdata")
#get threshold
sig.thr <- get.sig.thr(perms.1000.qtl.GFR.C2.192.Hi , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE )

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.GFR.C2.sex.Hi.perm1000.png", width = 2000, height = 1000, res = 100)
plot(Hi, sig.thr = sig.thr, sig.col = c("red", "orange", "chartreuse"), main = "Hi: Col4a5xDO logC2GFR Interaction only QTL map w/ 1000 perms")
dev.off()










