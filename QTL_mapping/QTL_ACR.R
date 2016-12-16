#QTL analsyses of ACR as covariate
#Two different analysis
#1 ACR with sex as covariate at 6wks, 10wks, and 15wks
#2 Albumin with creatinine as covariate at 6wks, 10wks, and 15wks

library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load sample
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
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
#head(pheno$Alb15WK_log[order(pheno$Alb15WK_log, decreasing=F)])

#check out distribution of ACR
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/hist.ACR.compare.png", width = 1000, height = 1000, res = 100)
layout(matrix(1:6,2,3))
hist(pheno$ACR6WK)
hist(pheno$ACR6WK_log)
hist(pheno$ACR10WK)
hist(pheno$ACR10WK_log)
hist(pheno$ACR15WK)
hist(pheno$ACR15WK_log)
dev.off()
#check distribution of Albumin 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/hist.Albumin.compare.png", width = 1000, height = 1000, res = 100)
layout(matrix(1:6,2,3))
hist(pheno$Alb6WK)
hist(pheno$Alb6WK_log)
hist(pheno$Alb10WK)
hist(pheno$Alb10WK_log)
hist(pheno$Alb15WK)
hist(pheno$Alb15WK_log)
dev.off()
#check distribution of Creatinine
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/hist.Creatinine.compare.png", width = 1000, height = 1000, res = 100)
layout(matrix(1:6,2,3))
hist(pheno$Creat6WK)
hist(pheno$Creat6WK_log)
hist(pheno$Creat10WK)
hist(pheno$Creat10WK_log)
hist(pheno$Creat15WK)
hist(pheno$Creat15WK_log)
dev.off()

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#create a subset separating Male and male samples
best.genoprobs.192 <- best.genoprobs.192[ order( rownames( best.genoprobs.192)),,]
pheno <- pheno[ order( rownames(pheno)),]

pheno.F <- subset(pheno, Sex == "F")
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
#Female 6wk
temp <- sex.covar.F
temp$creat6wk <- pheno.F$Creat6WK_log
covar.sex.creat6wk.F <- temp
#Female 10wk
temp <- sex.covar.F
temp$creat10wk <- pheno.F$Creat10WK_log
covar.sex.creat10wk.F <- temp
#Female 15wk
temp <- sex.covar.F
temp$creat15wk <- pheno.F$Creat15WK_log
covar.sex.creat15wk.F <- temp
#Male 6wk
temp <- sex.covar.M
temp$creat6wk <- pheno.M$Creat6WK_log
covar.sex.creat6wk.M <- temp
#Male 10wk
temp <- sex.covar.M
temp$creat10wk <- pheno.M$Creat10WK_log
covar.sex.creat10wk.M <- temp
#Male 15wk
temp <- sex.covar.M
temp$creat15wk <- pheno.M$Creat15WK_log
covar.sex.creat15wk.M <- temp

#kinship mapping
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
K_GS.F <- kinship.probs(genoprob.F, snps = GM_snps, bychr = TRUE)
K_GS.M <- kinship.probs(genoprob.M, snps = GM_snps, bychr = TRUE)

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.F.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.M.Rdata")

#QTL mapping of ACR with sex as covariate at 3 different time points
#Both sexes
#log ACR6WKS
qtl.log.ACR6WK.192 <- scanone( pheno = pheno, pheno.col = "ACR6WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.log.ACR6WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.Rdata")
#log ACR10WKS
qtl.log.ACR10WK.192 <- scanone( pheno = pheno, pheno.col = "ACR10WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.log.ACR10WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.Rdata")
#log ACR15WKS
qtl.log.ACR15WK.192 <- scanone( pheno = pheno, pheno.col = "ACR15WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.log.ACR15WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.Rdata")
#Females
#log ACR6WKS
qtl.log.ACR6WK.192.F <- scanone( pheno = pheno.F, pheno.col = "ACR6WK_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
save(qtl.log.ACR6WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.F.Rdata")
#log ACR10WKS
qtl.log.ACR10WK.192.F <- scanone( pheno = pheno, pheno.col = "ACR10WK_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
save(qtl.log.ACR10WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.F.Rdata")
#log ACR15WKS
qtl.log.ACR15WK.192.F <- scanone( pheno = pheno, pheno.col = "ACR15WK_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
save(qtl.log.ACR15WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.F.Rdata")
#Males
#log ACR6WKS
qtl.log.ACR6WK.192.M <- scanone( pheno = pheno.M, pheno.col = "ACR6WK_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
save(qtl.log.ACR6WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.M.Rdata")
#log ACR10WKS
qtl.log.ACR10WK.192.M <- scanone( pheno = pheno, pheno.col = "ACR10WK_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
save(qtl.log.ACR10WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.M.Rdata")
#log ACR15WKS
qtl.log.ACR15WK.192.M <- scanone( pheno = pheno, pheno.col = "ACR15WK_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
save(qtl.log.ACR15WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.M.Rdata")
#QTL mapping of Albumin with creatinine as covariate at 3 different time points
#log Albumin 6WKS
qtl.log.Alb6WK.192 <- scanone( pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat6wk, snps = GM_snps)
save(qtl.log.Alb6WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
#log Albumin 10WKS
qtl.log.Alb10WK.192 <- scanone( pheno = pheno, pheno.col = "Alb10WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat10wk, snps = GM_snps)
save(qtl.log.Alb10WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")
#log Albumin 15WKS
qtl.log.Alb15WK.192 <- scanone( pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat15wk, snps = GM_snps)
save(qtl.log.Alb15WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
#Females
#log Albumin 6WKS
qtl.log.Alb6WK.192.F <- scanone( pheno = pheno.F, pheno.col = "Alb6WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat6wk.F, snps = GM_snps)
save(qtl.log.Alb6WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.F.Rdata")
#log Albumin 10WKS
qtl.log.Alb10WK.192.F <- scanone( pheno = pheno.F, pheno.col = "Alb10WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat10wk.F, snps = GM_snps)
save(qtl.log.Alb10WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.F.Rdata")
#log Albumin 15WKS
qtl.log.Alb15WK.192.F <- scanone( pheno = pheno.F, pheno.col = "Alb15WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat15wk.F, snps = GM_snps)
save(qtl.log.Alb15WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.F.Rdata")
#Males
#log Albumin 6WKS
qtl.log.Alb6WK.192.M <- scanone( pheno = pheno.M, pheno.col = "Alb6WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat6wk.M, snps = GM_snps)
save(qtl.log.Alb6WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.M.Rdata")
#log Albumin 10WKS
qtl.log.Alb10WK.192.M <- scanone( pheno = pheno.M, pheno.col = "Alb10WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat10wk.M, snps = GM_snps)
save(qtl.log.Alb10WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.M.Rdata")
#log Albumin 15WKS
qtl.log.Alb15WK.192.M <- scanone( pheno = pheno.M, pheno.col = "Alb15WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat15wk.M, snps = GM_snps)
save(qtl.log.Alb15WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.M.Rdata")

#Run permutations 
#Perms 100: ACR
#log ACR6WKS
perms.100.qtl.log.ACR6WK.192  <- scanone.perm( pheno = pheno, pheno.col = "ACR6WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 100)
save(perms.100.qtl.log.ACR6WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.log.ACR6WK.192.Rdata" )
#log ACR10WKS
perms.100.qtl.log.ACR10WK.192  <- scanone.perm( pheno = pheno, pheno.col = "ACR10WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 100)
save(perms.100.qtl.log.ACR10WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.log.ACR10WK.192.Rdata" )
#log ACR15WKS
perms.100.qtl.log.ACR15WK.192  <- scanone.perm( pheno = pheno, pheno.col = "ACR15WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 100)
save(perms.100.qtl.log.ACR15WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm100/perms.100.qtl.log.ACR15WK.192.Rdata" )
#Perms 1000: ACR
#log ACR6WKS
perms.1000.qtl.log.ACR6WK.192  <- scanone.perm( pheno = pheno, pheno.col = "ACR6WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.ACR6WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.ACR6WK.192.Rdata" )
#log ACR10WKS
perms.1000.qtl.log.ACR10WK.192  <- scanone.perm( pheno = pheno, pheno.col = "ACR10WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.ACR10WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.ACR10WK.192.Rdata" )
#log ACR15WKS
perms.1000.qtl.log.ACR15WK.192  <- scanone.perm( pheno = pheno, pheno.col = "ACR15WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.ACR15WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.ACR15WK.192.Rdata" )
#Perms 1000: Albumin with creat covar
#log Alb6WKS
perms.1000.qtl.log.Alb6WK.192  <- scanone.perm( pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, addcovar = covar.sex.creat6wk, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.Alb6WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.192.Rdata" )
#log Alb10WKS
perms.1000.qtl.log.Alb10WK.192  <- scanone.perm( pheno = pheno, pheno.col = "Alb10WK_log", probs = best.genoprobs.192, addcovar = covar.sex.creat10wk, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.Alb10WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.192.Rdata" )
#log Alb15WKS
perms.1000.qtl.log.Alb15WK.192  <- scanone.perm( pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, addcovar = covar.sex.creat15wk, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.Alb15WK.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15WK.192.Rdata" )

#Assign threshold
#threshold: ACR 1000
thr.1000.qtl.logACR.6wk <- get.sig.thr( perms.1000.qtl.log.ACR6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logACR.10wk <- get.sig.thr( perms.1000.qtl.log.ACR10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logACR.15wk <- get.sig.thr( perms.1000.qtl.log.ACR15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#threshold: Alb 1000
thr.1000.qtl.logAlb.6wk <- get.sig.thr( perms.1000.qtl.log.Alb6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.10wk <- get.sig.thr( perms.1000.qtl.log.Alb10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.15wk <- get.sig.thr( perms.1000.qtl.log.Alb15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)


#Create plots
#plot QTL: ACR 
#ACR : both sexes
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perm1000.log.ACR6wk.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR6WK.192, sig.thr = thr.1000.qtl.logACR.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.ACR.6WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perm1000.log.ACR10wk.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR10WK.192, sig.thr = thr.1000.qtl.logACR.10wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.ACR.10WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perm1000.log.ACR15wk.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR15WK.192, sig.thr = thr.1000.qtl.logACR.15wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.ACR.15WK QTL perms.1000")
dev.off()
#ACR : Females no threshold
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.ACR6wk.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR6WK.192.F, main = "Col4a5xDO log.ACR.6WK QTL Females")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.ACR10wk.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR10WK.192.F, main = "Col4a5xDO log.ACR.10WK QTL Females")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.ACR15wk.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR15WK.192.F, main = "Col4a5xDO log.ACR.15WK QTL Females")
dev.off()
#ACR : Males no threshold 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.ACR6wk.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR6WK.192.M, main = "Col4a5xDO log.ACR.6WK QTL Males")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.ACR10wk.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR10WK.192.M, main = "Col4a5xDO log.ACR.10WK QTL Males")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.ACR15wk.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.ACR15WK.192.M, main = "Col4a5xDO log.ACR.15WK QTL Males")
dev.off()
#Create plots
#plot QTL: Alb
#Alb : Both sexes
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perm1000.log.Alb6wk.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192, sig.thr = thr.1000.qtl.logAlb.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.6WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perm1000.log.Alb10wk.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192, sig.thr = thr.1000.qtl.logAlb.10wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.10WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perm1000.log.Alb15wk.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192, sig.thr = thr.1000.qtl.logAlb.15wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.15WK QTL perms.1000")
dev.off()
#Alb : Females
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb6wk.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192.F, main = "Col4a5xDO log.Alb.6WK QTL Females")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb10wk.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192.F, main = "Col4a5xDO log.Alb.10WK QTL Females")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb15wk.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192.F, main = "Col4a5xDO log.Alb.15WK QTL Females")
dev.off()
#Alb : Males
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb6wk.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192.M, main = "Col4a5xDO log.Alb.6WK QTL Males")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb10wk.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192.M, main = "Col4a5xDO log.Alb.10WK QTL Males")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb15wk.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192.M, main = "Col4a5xDO log.Alb.15WK QTL Males")
dev.off()





