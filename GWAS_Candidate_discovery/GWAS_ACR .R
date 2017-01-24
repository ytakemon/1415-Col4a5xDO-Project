library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load sample
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR10WK_log <- log(pheno$ACR10WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)

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

#create a subset separating Male and male samples
best.genoprobs.192 <- best.genoprobs.192[ order( rownames( best.genoprobs.192)),,]
pheno <- pheno[ order( rownames(pheno)),]

pheno.F <- subset(pheno, Sex =="F")
pheno.M <- subset(pheno, Sex == "M")

genoprob.F <- best.genoprobs.192[rownames(pheno.F),,]
genoprob.M <- best.genoprobs.192[rownames(pheno.M),,]

sex.covar$MouseID <- rownames(sex.covar)
sex.covar <- sex.covar[ order( sex.covar$MouseID),]

sex.covar.F <- sex.covar[rownames(pheno.F),]
sex.covar.M <- sex.covar[rownames(pheno.M),]

sex.covar$MouseID <- NULL
sex.covar.F$MouseID <- NULL
sex.covar.M$MouseID <- NULL

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

#kinship mapping
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.F.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.M.Rdata")

#creating manual permuations
#############################
#males <- which(sex.covar[,1] == 1)
#females <- which(sex.covar[,1] == 0)

#males.perm <- sample( rownames( pheno)[males])
#females.perm <- sample( rownames(pheno)[females])

#rownames(pheno)[males] <- males.perm 
#rownames(pheno)[females] <- females.perm 

#rownames(sex.covar)[males] <- males.perm
#rownames(sex.covar)[females] <- females.perm
#############################

sex.covar <- as.matrix(sex.covar)
covar.sex.creat6wk <- as.matrix(covar.sex.creat6wk)
covar.sex.creat10wk <- as.matrix(covar.sex.creat10wk)
covar.sex.creat15wk <- as.matrix(covar.sex.creat15wk)
#GWAS ACR 6
GWAS.ACR6WK <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = sex.covar, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.ACR6WK, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.ACR6wk.Rdata")

#GWAS ACR 10
GWAS.ACR10WK <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "ACR10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = sex.covar, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.ACR10WK, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.ACR10wk.Rdata")

#GWAS ACR 15
GWAS.ACR15WK <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "ACR15WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = sex.covar, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.ACR15WK, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.ACR15wk.Rdata")

#GWAS ALB 6
GWAS.Alb6WK <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = covar.sex.creat6wk, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.Alb6WK, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb6WK.Rdata")

#GWAS ALB 10
GWAS.Alb10WK <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "Alb10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = covar.sex.creat10wk, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.Alb10WK, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb10wk.Rdata")

#GWAS ALB 15
pheno[pheno ==  -Inf] = NA
GWAS.Alb15WK <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "Alb15WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = covar.sex.creat15wk, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.Alb15WK, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb15wk.Rdata")

#create permutation then run this 1000 times
#ACR6
min.p <- min(sapply(GWAS.ACR6WK[1:19], function(z){min(z$p.value)}))
write( min.p, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/ACR6.perm/perm.",args[1],".txt"))
#ACR10
min.p <- min(sapply(GWAS.ACR10WK[1:19], function(z){min(z$p.value)}))
write( min.p, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/ACR10.perm/perm.",args[1],".txt"))
#ACR15
min.p <- min(sapply(GWAS.ACR15WK[1:19], function(z){min(z$p.value)}))
write( min.p, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/ACR15.perm/perm.",args[1],".txt"))
#Alb6
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb6WK.Rdata")
min.p <- min(sapply(GWAS.Alb6WK[1:19], function(z){min(z$p.value)}))
write( min.p, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/Alb6.perm/perm.",args[1],".txt"))
#Alb10
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb10wk.Rdata")
min.p <- min(sapply(GWAS.Alb10WK[1:19], function(z){min(z$p.value)}))
write( min.p, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/Alb10.perm/perm.",args[1],".txt"))
#Alb15
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb15wk.Rdata")
min.p <- min(sapply(GWAS.Alb15WK[1:19], function(z){min(z$p.value)}))
write( min.p, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/Alb15.perm/perm.",args[1],".txt"))

#jump to cadillac folder with all perm files
grep -h . *.txt > all.perm.txt

#read into R
ACR6.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/ACR6.perm/all.perm.txt", sep ="\t", header = F )
ACR10.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/ACR10.perm/all.perm.txt", sep ="\t", header = F )
ACR15.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/ACR15.perm/all.perm.txt", sep ="\t", header = F )
Alb6.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/Alb6.perm/all.perm.txt", sep ="\t", header = F )
Alb10.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/Alb10.perm/all.perm.txt", sep ="\t", header = F )
Alb15.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/Alb15.perm/all.perm.txt", sep ="\t", header = F )
#load gwas file
load( "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.ACR6wk.Rdata")
load( "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.ACR10wk.Rdata")
load( "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.ACR15wk.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb6WK.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb10wk.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.Alb15wk.Rdata")

#thresholding
gwas.ACR6.perm <- get.sig.thr( -log10(ACR6.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.ACR10.perm <- get.sig.thr( -log10(ACR10.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.ACR15.perm <- get.sig.thr( -log10(ACR15.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.Alb6.perm <- get.sig.thr( -log10(Alb6.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.Alb10.perm <- get.sig.thr( -log10(Alb10.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.Alb15.perm <- get.sig.thr( -log10(Alb15.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

#plots
#ACR6
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.ACR6wk.png", width = 2000, height = 1000, res = 100)
plot(GWAS.ACR6WK, sig.thr = gwas.ACR6.perm, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of ACR6WK")
dev.off()
#ACR10
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.ACR10wk.png", width = 2000, height = 1000, res = 100)
plot(GWAS.ACR10WK, ylim = c(0 ,max(gwas.ACR10.perm)), sig.thr = gwas.ACR10.perm, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of ACR10WK")
dev.off()
#ACR15
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.ACR15wk.png", width = 2000, height = 1000, res = 100)
plot(GWAS.ACR15WK, ylim = c(0 ,max(gwas.ACR15.perm)), sig.thr = gwas.ACR15.perm, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of ACR15WK")
dev.off()
#Alb6
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.Alb6wk.png", width = 2000, height = 1000, res = 100)
plot(GWAS.Alb6WK, sig.thr = gwas.Alb6.perm, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of Alb6WK")
dev.off()
#Alb10
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.Alb10wk.png", width = 2000, height = 1000, res = 100)
plot(GWAS.Alb10WK, ylim = c(0 ,max(gwas.Alb10.perm)),sig.thr = gwas.Alb10.perm, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of Alb10WK")
dev.off()
#Alb15
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.Alb15wk.png", width = 2000, height = 1000, res = 100)
plot(GWAS.Alb15WK, ylim = c(0 ,max(gwas.Alb15.perm)), sig.thr = gwas.Alb15.perm, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of Alb15WK")
dev.off()

#################################################################
#Coefficient Allele effect plot

########################################6WK
#Load QTL data: 6WKS
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.Rdata")
ACR_qtl <- qtl.log.ACR6WK.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
Alb_qtl <- qtl.log.Alb6WK.192

##############Chr2
#ACR
ACR_qtl <- qtl.log.ACR6WK.192
ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR6.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 2, main = "Chr 2 ACR6WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb6WK.192
Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb6.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 2, main = "Chr 2 Alb6WK_log")
dev.off()

##############Chr11
#ACR
ACR_qtl <- qtl.log.ACR6WK.192
#ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR6.chr11.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 11, main = "Chr 11 ACR6WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb6WK.192
#Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb6.chr11.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 11, main = "Chr 11 Alb6WK_log")
dev.off()

##############Chr12
#ACR
ACR_qtl <- qtl.log.ACR6WK.192
#ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR6.chr12.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 12, main = "Chr 12 ACR6WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb6WK.192
#Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb6.chr12.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 12, main = "Chr 12 Alb6WK_log")
dev.off()


########################################10WK
#Load QTL data: 10WKS
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.Rdata")
ACR_qtl <- qtl.log.ACR10WK.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")
Alb_qtl <- qtl.log.Alb10WK.192

##############Chr2
#ACR
ACR_qtl <- qtl.log.ACR10WK.192
ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR10.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 2, main = "Chr 2 ACR10WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb10WK.192
Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 2, main = "Chr 2 Alb10WK_log")
dev.off()

##############Chr4
#ACR
ACR_qtl <- qtl.log.ACR10WK.192
#ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR10.chr4.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 4, main = "Chr 4 ACR10WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb10WK.192
#Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr4.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 4, main = "Chr 4 Alb10WK_log")
dev.off()

##############Chr11
#ACR
ACR_qtl <- qtl.log.ACR10WK.192
#ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR10.chr11.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 11, main = "Chr 11 ACR10WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb10WK.192
#Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr11.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 11, main = "Chr 11 Alb10WK_log")
dev.off()

##############Chr13
#ACR
ACR_qtl <- qtl.log.ACR10WK.192
ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR10.chr13.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 13, main = "Chr 13 ACR10WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb10WK.192
Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr13.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 13, main = "Chr 13 Alb10WK_log")
dev.off()

########################################15WK
#Load QTL data: 15WKS
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.Rdata")
ACR_qtl <- qtl.log.ACR15WK.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
Alb_qtl <- qtl.log.Alb15WK.192

##############Chr15
#ACR
ACR_qtl <- qtl.log.ACR15WK.192
#ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR15.chr15.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 15, main = "Chr 15 ACR15WK_log")
dev.off()
#Alb
Alb_qtl <- qtl.log.Alb15WK.192
#Alb_qtl$coef$A[abs(Alb_qtl$coef$A) > 1.5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb15.chr15.png", width = 1500, height = 1000, res = 100)
coefplot(Alb_qtl , chr = 15, main = "Chr 15 Alb15WK_log")
dev.off()


######################################################################################################################################################
######################################################################################################################################################
####################################################################CANDIDATE GENES###################################################################

####################################################################6 WEEKS###########################################################################
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.Rdata")
ACR_qtl <- qtl.log.ACR6WK.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
Alb_qtl <- qtl.log.Alb6WK.192
#####################################################################6WK: Chr2
#ACR6: CHR2
interval = bayesint(ACR_qtl, chr = 2)
knitr::kable(interval)
#|            |marker      |chr |        pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|----------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|UNCHS003833 |UNCHS003833 |2   |   3.176721|  0.000312|  9.702188| 16.73734| 3.634467| 0.0191707|    1.717363|
#|UNCHS007859 |UNCHS007859 |2   | 173.869006| 87.336158| 14.486049| 25.66447| 5.572968| 0.0005781|    3.237980|
#|JAX00512278 |JAX00512278 |2   | 174.374458| 87.751661|  9.171967| 15.77716| 3.425968| 0.0272321|    1.564920|
#pt1
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 2, 
	start = 110,
	end = 114, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr2.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt2
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 2, 
	start = 172,
	end = 176, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr2.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 3, show.sdps = TRUE)
dev.off()
#Alb6: CHR2
#pt1
interval = bayesint(Alb_qtl, chr = 2)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|---------:|--------:|--------:|---------:|-----------:|
#|JAX00499427 |JAX00499427 |2   | 109.1249| 49.54448| 10.943882| 19.23998| 4.177908| 0.0074679|    2.126801|
#|JAX00098823 |JAX00098823 |2   | 113.4831| 50.68170| 16.654300| 30.24075| 6.566695| 0.0000858|    4.066649|
#|UNC4549601  |UNC4549601  |2   | 174.0296| 87.67960|  9.869729| 17.24974| 3.745734| 0.0158549|    1.799837|
#pt1
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 2, 
	start = 110,
	end = 114, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr2.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 5, show.sdps = TRUE)
dev.off()
#pt2
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 2, 
	start = 172,
	end = 176, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr2.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 3, show.sdps = TRUE)
dev.off()

#####################################################################6WK: Chr11
#ACR6: CHR11
interval = bayesint(ACR_qtl, chr = 11)
knitr::kable(interval)
#|            |marker      |chr |       pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|JAX00303276 |JAX00303276 |11  |  6.953372|  2.298799| 10.506113| 18.20398| 3.952944| 0.0110828|    1.955350|
#|JAX00317684 |JAX00317684 |11  | 89.411925| 48.159751| 16.775678| 30.11541| 6.539478| 0.0000904|    4.043636|
#|UNC20099016 |UNC20099016 |11  | 90.033534| 49.354168|  9.888831| 17.07668| 3.708153| 0.0169084|    1.771898|
#pt1
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 11, 
	start = 4,
	end = 11, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr11.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt2
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 11, 
	start = 70,
	end = 76, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr11.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt3
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 11, 
	start = 87,
	end = 90, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr11.pt3.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 5, show.sdps = TRUE)
dev.off()
#Alb6: CHR11
interval = bayesint(Alb_qtl, chr = 11)
knitr::kable(interval)
#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|JAX00303381 |JAX00303381 |11  |  7.359184|  2.34442| 12.44737| 22.06639| 4.791655| 0.0024743|    2.606544|
#|JAX00317690 |JAX00317690 |11  | 89.429042| 48.19158| 17.95163| 32.84497| 7.132195| 0.0000283|    4.548272|
#|UNC20099016 |UNC20099016 |11  | 90.033534| 49.35417| 10.31149| 18.06537| 3.922846| 0.0116785|    1.932613|
#pt1
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 11, 
	start = 4,
	end = 11, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr11.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt2
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 11, 
	start = 70,
	end = 76, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr11.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt3
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 11, 
	start = 87,
	end = 90, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr11.pt3.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 6, show.sdps = TRUE)
dev.off()

#####################################################################6WK: Chr12
#ACR6: CHR12
interval = bayesint(ACR_qtl, chr = 12)
knitr::kable(interval)
#|            |marker      |chr |        pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|----------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|UNC20604636 |UNC20604636 |12  |   9.471158|  2.282116|  8.584962| 14.72067| 3.196553| 0.0397514|    1.400647|
#|UNC21410263 |UNC21410263 |12  |  74.675566| 29.685522| 14.163369| 25.04679| 5.438842| 0.0007444|    3.128169|
#|UNC21973760 |UNC21973760 |12  | 116.288534| 59.935034|  8.131618| 13.90938| 3.020382| 0.0528169|    1.277227|
#pt1
chr <- 12
chr12.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 12, 
	start = 70,
	end = 82, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr12.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr12.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt2
chr <- 12
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 12, 
	start = 107,
	end = 109.5, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6.Chr12.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr12.genes, thr = 3, show.sdps = TRUE)
dev.off()
#Alb6: CHR12
interval = bayesint(Alb_qtl, chr = 12)
knitr::kable(interval)
#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|JAX00337821 |JAX00337821 |12  |  72.95626| 28.52262| 10.02364| 17.53345| 3.807340| 0.0142620|    1.845819|
#|UNCHS033908 |UNCHS033908 |12  |  74.01750| 29.60898| 14.74737| 26.48551| 5.751255| 0.0004123|    3.384761|
#|UNCHS034675 |UNCHS034675 |12  | 109.09296| 54.61919| 10.29093| 18.02731| 3.914581| 0.0118474|    1.926378|
#pt1
chr <- 12
chr12.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 12, 
	start = 70,
	end = 82, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr12.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr12.genes, thr = 3, show.sdps = TRUE)
dev.off()
#pt2
chr <- 12
chr12.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb6WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk, 
	snps = GM_snps, 
	chr = 12, 
	start = 107,
	end = 109.5, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr12.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr12.genes, thr = 3, show.sdps = TRUE)
dev.off()


####################################################################10 WEEKS###########################################################################
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.Rdata")
ACR_qtl <- qtl.log.ACR10WK.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")
Alb_qtl <- qtl.log.Alb10WK.192
#####################################################################10WK: Chr2
#ACR10: CHR2
interval = bayesint(ACR_qtl, chr = 2)
knitr::kable(interval)
#|            |marker      |chr |        pos|        cM|  perc.var|       lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|----------:|---------:|---------:|---------:|--------:|---------:|-----------:|
#|UNCHS003833 |UNCHS003833 |2   |   3.176721|  0.000312|  7.502436| 12.478061| 2.709576| 0.0858936|   1.0660393|
#|JAX00506307 |JAX00506307 |2   | 145.703227| 63.905578| 10.401685| 17.573387| 3.816012| 0.0140505|   1.8523094|
#|UNC4609888  |UNC4609888  |2   | 181.998284| 93.369576|  4.586822|  7.512557| 1.631331| 0.3775282|   0.4230506|
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 2, 
	start = 101,
	end = 114, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10.Chr2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 4, show.sdps = TRUE)
dev.off()
#Alb10: CHR2
interval = bayesint(Alb_qtl, chr = 2)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNCHS005963 |UNCHS005963 |2   | 101.7360| 45.94599| 11.89948| 20.39737| 4.429233| 0.0047724|    2.321265|
#|JAX00098817 |JAX00098817 |2   | 113.4217| 50.67787| 16.48704| 29.00710| 6.298813| 0.0001443|    3.840867|
#|UNC3776734  |UNC3776734  |2   | 113.6730| 50.81143| 11.95717| 20.50284| 4.452134| 0.0045801|    2.339128|
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat10wk, 
	snps = GM_snps, 
	chr = 2, 
	start = 101,
	end = 114, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb10.Chr2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 4, show.sdps = TRUE)
dev.off()

#####################################################################10WK: Chr4
#ACR10: CHR4
interval = bayesint(ACR_qtl, chr = 4)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNC7389537  |UNC7389537  |4   |  64.8536| 28.64806| 12.32392| 21.04337| 4.569509| 0.0037065|    2.431037|
#|UNCHS012670 |UNCHS012670 |4   | 122.8761| 49.38360| 15.24768| 26.46993| 5.747872| 0.0004150|    3.381967|
#|UNC8168429  |UNC8168429  |4   | 126.6525| 51.18213| 11.15567| 18.92551| 4.109621| 0.0084239|    2.074485|
#pt1
chr <- 4
chr4.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 4, 
	start = 60,
	end = 75, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10.Chr4.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr4.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt2
chr <- 4
chr4.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 4, 
	start = 118,
	end = 126, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10.Chr4.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr4.genes, thr = 4, show.sdps = TRUE)
dev.off()
#Alb10: CHR4
interval = bayesint(Alb_qtl, chr = 4)
knitr::kable(interval)
#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNCHS011502 |UNCHS011502 |4   |  63.27366| 27.77400| 12.62087| 21.72110| 4.716677| 0.0028374|    2.547087|
#|JAX00124064 |JAX00124064 |4   | 122.92428| 49.40889| 15.97110| 28.01551| 6.083491| 0.0002185|    3.660600|
#|JAX00564645 |JAX00564645 |4   | 124.55930| 49.65606| 13.30286| 22.98263| 4.990615| 0.0017165|    2.765351|
#pt1
chr <- 4
chr4.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat10wk, 
	snps = GM_snps, 
	chr = 4, 
	start = 60,
	end = 75, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb10.Chr4.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr4.genes, thr = 4, show.sdps = TRUE)
dev.off()
#pt2
chr <- 4
chr4.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat10wk, 
	snps = GM_snps, 
	chr = 4, 
	start = 118,
	end = 126, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb10.Chr4.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr4.genes, thr = 4, show.sdps = TRUE)
dev.off()

#####################################################################10WK: Chr11
#ACR10: CHR11
interval = bayesint(ACR_qtl, chr = 11)
knitr::kable(interval)
#|            |marker      |chr |        pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|----------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|JAX00302685 |JAX00302685 |11  |   3.263868|  0.020216|  7.134328| 11.84258| 2.571584| 0.1058560|   0.9752845|
#|UNCHS029750 |UNCHS029750 |11  |   7.504426|  2.372091| 13.609066| 23.40599| 5.082546| 0.0014481|   2.8392143|
#|UNCHS032292 |UNCHS032292 |11  | 109.993626| 63.634193|  6.263726| 10.34959| 2.247385| 0.1696128|   0.7705415|
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 11, 
	start = 6,
	end = 10, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10.Chr11.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 4, show.sdps = TRUE)
dev.off()
#Alb10: CHR11
interval = bayesint(Alb_qtl, chr = 11)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|---------:|--------:|--------:|---------:|-----------:|
#|UNC19018757 |UNC19018757 |11  | 6.926399| 2.295767|  9.601176| 16.25117| 3.528896| 0.0229182|    1.639820|
#|UNCHS029750 |UNCHS029750 |11  | 7.504426| 2.372091| 15.668772| 27.43729| 5.957932| 0.0002780|    3.556014|
#|UNC19055538 |UNC19055538 |11  | 9.394088| 3.605907| 10.505942| 17.87067| 3.880567| 0.0125673|    1.900757|
chr <- 11
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat10wk, 
	snps = GM_snps, 
	chr = 11, 
	start = 6,
	end = 10, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb10.Chr11.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 4, show.sdps = TRUE)
dev.off()

#####################################################################10WK: Chr13
#ACR10: CHR13
interval = bayesint(ACR_qtl, chr = 13)
knitr::kable(interval)
#|            |marker      |chr |       pos|       cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|---------:|--------:|--------:|---------:|-----------:|
#|UNC22439783 |UNC22439783 |13  |  37.86027| 12.91926|  8.276310| 13.82232| 3.001479| 0.0544333|    1.264135|
#|UNC23234166 |UNC23234166 |13  | 101.78100| 49.09407| 14.953252| 25.91506| 5.627383| 0.0005216|    3.282683|
#|UNCHS037276 |UNCHS037276 |13  | 118.44275| 62.66769|  8.494274| 14.20298| 3.084138| 0.0476867|    1.321603|
chr <- 13
chr13.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 13, 
	start = 100,
	end = 105, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10.Chr13.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr13.genes, thr = 3, show.sdps = TRUE)
dev.off()
#Alb10: CHR13
interval = bayesint(Alb_qtl, chr = 13)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|---------:|--------:|--------:|---------:|-----------:|
#|UNC19018757 |UNC19018757 |11  | 6.926399| 2.295767|  9.601176| 16.25117| 3.528896| 0.0229182|    1.639820|
#|UNCHS029750 |UNCHS029750 |11  | 7.504426| 2.372091| 15.668772| 27.43729| 5.957932| 0.0002780|    3.556014|
#|UNC19055538 |UNC19055538 |11  | 9.394088| 3.605907| 10.505942| 17.87067| 3.880567| 0.0125673|    1.900757|
chr <- 13
chr11.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb10WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat10wk, 
	snps = GM_snps, 
	chr = 13, 
	start = 100,
	end = 105, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb10.Chr13.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr11.genes, thr = 3, show.sdps = TRUE)
dev.off()


####################################################################15 WEEKS###########################################################################
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.Rdata")
ACR_qtl <- qtl.log.ACR15WK.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
Alb_qtl <- qtl.log.Alb15WK.192
#####################################################################10WK: Chr2
#ACR15: CHR15
interval = bayesint(ACR_qtl, chr = 15)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|       p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|-------:|-----------:|
#|UNC26049729 |UNC26049729 |15  | 87.49891| 37.34539| 21.05914| 39.96363| 8.677991| 1.3e-06|    5.893086|
#|UNCHS041222 |UNCHS041222 |15  | 87.98976| 37.74115| 22.72292| 43.56361| 9.459718| 3.0e-07|    6.585735|
#|JAX00408925 |JAX00408925 |15  | 88.23271| 38.10487| 20.99527| 39.82695| 8.648313| 1.4e-06|    5.866937|
chr <- 15
chr15.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="ACR15WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 15, 
	start = 85,
	end = 90, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR15.Chr15.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr15.genes, thr = 2, show.sdps = TRUE)
dev.off()
#Alb15: CHR15
interval = bayesint(Alb_qtl, chr = 15)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|       lod|     p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|---------:|-----:|-----------:|
#|UNC26048887 |UNC26048887 |15  | 87.44608| 37.22123| 22.18870| 42.39931|  9.206894| 4e-07|    6.360931|
#|UNCHS041212 |UNCHS041212 |15  | 87.67108| 37.35712| 24.60839| 47.73815| 10.366207| 0e+00|    7.397302|
#|UNC26057169 |UNC26057169 |15  | 88.00011| 37.75250| 22.22703| 42.48259|  9.224977| 4e-07|    6.376986|
#Chr15

chr <- 15
chr15.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="Alb15WK_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = covar.sex.creat15wk, 
	snps = GM_snps, 
	chr = 15, 
	start = 85,
	end = 90, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb15.Chr15.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr15.genes, thr = 3, show.sdps = TRUE)
dev.off()


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
##############################################################Check Genotype Distribution##################################################################
#Alb 6wks
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
qtl <- qtl.log.Alb6WK.192
#Chr2
interval = bayesint(qtl, chr = 2)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|---------:|--------:|--------:|---------:|-----------:|
#|JAX00499427 |JAX00499427 |2   | 109.1249| 49.54448| 10.943882| 19.23998| 4.177908| 0.0074679|    2.126801|
#|JAX00098823 |JAX00098823 |2   | 113.4831| 50.68170| 16.654300| 30.24075| 6.566695| 0.0000858|    4.066649|
#|UNC4549601  |UNC4549601  |2   | 174.0296| 87.67960|  9.869729| 17.24974| 3.745734| 0.0158549|    1.799837|
GM_2 <- GM_snps[GM_snps$chr == 2,]
GM_2.113 <- GM_2[GM_2$pos > 113,]
GM_2.114 <- GM_2.113[GM_2.113$pos < 114,]
#UNC3769267		113.0337
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb6.chr2.geno.dist.113.0.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, snp.id = "UNC3769267", snps = GM_snps)
dev.off()
#JAX00098823	113.4831
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb6.chr2.geno.dist.113.4.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, snp.id = "JAX00098823", snps = GM_snps)
dev.off()
#UNC3779021		113.9611
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb6.chr2.geno.dist.113.9.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, snp.id = "UNC3779021", snps = GM_snps)
dev.off()
#Chr11
interval = bayesint(qtl, chr = 11)
knitr::kable(interval)
#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|JAX00303381 |JAX00303381 |11  |  7.359184|  2.34442| 12.44737| 22.06639| 4.791655| 0.0024743|    2.606544|
#|JAX00317690 |JAX00317690 |11  | 89.429042| 48.19158| 17.95163| 32.84497| 7.132195| 0.0000283|    4.548272|
#|UNC20099016 |UNC20099016 |11  | 90.033534| 49.35417| 10.31149| 18.06537| 3.922846| 0.0116785|    1.932613|
GM_11 <- GM_snps[GM_snps$chr == 11,]
GM_11.89.5 <- GM_11[GM_11$pos < 89.5,]
GM_11.89.3 <- GM_11.89.5[GM_11.89.5$pos > 89.3,]
#UNC20088171		89.30357
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb6.chr11.geno.dist.89.3.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, snp.id = "UNC20088171", snps = GM_snps)
dev.off()
#JAX00317690		89.429042
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb6.chr11.geno.dist.89.42.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, snp.id = "JAX00317690", snps = GM_snps)
dev.off()
#UNC20090864		89.49846
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb6.chr11.geno.dist.89.49.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, snp.id = "UNC20090864", snps = GM_snps)
dev.off()

#######################################################################################################################################################

#Alb 15wks
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
qtl <- qtl.log.Alb15WK.192
#Chr15
interval = bayesint(qtl, chr = 15)
knitr::kable(interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|       lod|     p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|---------:|-----:|-----------:|
#|UNC26048887 |UNC26048887 |15  | 87.44608| 37.22123| 22.18870| 42.39931|  9.206894| 4e-07|    6.360931|
#|UNCHS041212 |UNCHS041212 |15  | 87.67108| 37.35712| 24.60839| 47.73815| 10.366207| 0e+00|    7.397302|
#|UNC26057169 |UNC26057169 |15  | 88.00011| 37.75250| 22.22703| 42.48259|  9.224977| 4e-07|    6.376986|
GM_15 <- GM_snps[GM_snps$chr == 15,]
GM_15.89 <- GM_15[GM_15$pos < 88.2,]
GM_15.87 <- GM_15.89[GM_15.89$pos > 87.7,]
#Chck all three intervals

#UNC26048887 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb15.chr15.geno.dist.87.4.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, snp.id = "UNC26048887", snps = GM_snps)
dev.off()
#UNCHS041212
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb15.chr15.geno.dist.87.6.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, snp.id = "UNCHS041212", snps = GM_snps)
dev.off()
#UNC26057169
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/alb15.chr15.geno.dist.88.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, snp.id = "UNC26057169", snps = GM_snps)
dev.off()




