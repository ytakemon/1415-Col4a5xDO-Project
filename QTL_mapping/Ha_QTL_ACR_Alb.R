## QTL analysis for Ha model of ACR and Albumin
## Yuka Takemon
## Created 10/03/16
## Updataed 03/08/20

#	Two different analysis
#	ACR with sex as covariate at 6wks, 10wks, and 15wks
#	Albumin with sex and creatinine as covariate at 6wks, 10wks, and 15wks

library(DOQTL)
library(abind)
library(knitr)

dir <- "/projects/marralab/ytakemon_prj/Col4a5/"

#load sample
load(paste0(dir,"Data/consolidated/best.genoprobs.192.Rdata"))
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#read and clean up phenotype data
pheno <- read.delim(paste0(dir,"Data/consolidated/Phenotype/1415_master_pheno.txt"), sep = "\t", header = TRUE)
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
#K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
saveRDS(K_GS, paste0(dir,"Data/consolidated/K_GS.rds"))
#K_GS.F <- kinship.probs(genoprob.F, snps = GM_snps, bychr = TRUE)
#K_GS.M <- kinship.probs(genoprob.M, snps = GM_snps, bychr = TRUE)

K_GS <- readRDS(paste0(dir,"Data/consolidated/K_GS.rds"))

#QTL mapping of ACR with sex as covariate at 3 different time points
#Both sexes
#log ACR6WKS
qtl.log.ACR6WK.192 <- scanone( pheno = pheno, pheno.col = "ACR6WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.log.ACR6WK.192, file = paste0(dir,"Data/consolidated/qtl.log.ACR6WK.192.Rdata.rds"))
#log ACR10WKS
qtl.log.ACR10WK.192 <- scanone( pheno = pheno, pheno.col = "ACR10WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.log.ACR10WK.192, file = paste0(dir,"Data/consolidated/qtl.log.ACR10WK.192.Rdata.rds"))
#log ACR15WKS
qtl.log.ACR15WK.192 <- scanone( pheno = pheno, pheno.col = "ACR15WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.log.ACR15WK.192, file = paste0(dir,"Data/consolidated/qtl.log.ACR15WK.192.Rdata.rds"))

# #Females
# #log ACR6WKS
# qtl.log.ACR6WK.192.F <- scanone( pheno = pheno.F, pheno.col = "ACR6WK_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
# save(qtl.log.ACR6WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.F.Rdata")
# #log ACR10WKS
# qtl.log.ACR10WK.192.F <- scanone( pheno = pheno, pheno.col = "ACR10WK_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
# save(qtl.log.ACR10WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.F.Rdata")
# #log ACR15WKS
# qtl.log.ACR15WK.192.F <- scanone( pheno = pheno, pheno.col = "ACR15WK_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
# save(qtl.log.ACR15WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.F.Rdata")
# #Males
# #log ACR6WKS
# qtl.log.ACR6WK.192.M <- scanone( pheno = pheno.M, pheno.col = "ACR6WK_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
# save(qtl.log.ACR6WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.M.Rdata")
# #log ACR10WKS
# qtl.log.ACR10WK.192.M <- scanone( pheno = pheno, pheno.col = "ACR10WK_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
# save(qtl.log.ACR10WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.M.Rdata")
# #log ACR15WKS
# qtl.log.ACR15WK.192.M <- scanone( pheno = pheno, pheno.col = "ACR15WK_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
# save(qtl.log.ACR15WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR15WK.192.M.Rdata")

#QTL mapping of Albumin with creatinine as covariate at 3 different time points
#log Albumin 6WKS
qtl.log.Alb6WK.192 <- scanone( pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat6wk, snps = GM_snps)
save(qtl.log.Alb6WK.192, file = paste0(dir,"Data/consolidated/qtl.log.Alb6WK.192.Rdata.rds"))
#log Albumin 10WKS
qtl.log.Alb10WK.192 <- scanone( pheno = pheno, pheno.col = "Alb10WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat10wk, snps = GM_snps)
save(qtl.log.Alb10WK.192, file = paste0(dir,"Data/consolidated/qtl.log.Alb10WK.192.Rdata.rds"))
#log Albumin 15WKS
qtl.log.Alb15WK.192 <- scanone( pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat15wk, snps = GM_snps)
save(qtl.log.Alb15WK.192, file = paste0(dir,"Data/consolidated/qtl.log.Alb15WK.192.Rdata.rds"))

# #Females
# #log Albumin 6WKS
# qtl.log.Alb6WK.192.F <- scanone( pheno = pheno.F, pheno.col = "Alb6WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat6wk.F, snps = GM_snps)
# save(qtl.log.Alb6WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.F.Rdata")
# #log Albumin 10WKS
# qtl.log.Alb10WK.192.F <- scanone( pheno = pheno.F, pheno.col = "Alb10WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat10wk.F, snps = GM_snps)
# save(qtl.log.Alb10WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.F.Rdata")
# #log Albumin 15WKS
# qtl.log.Alb15WK.192.F <- scanone( pheno = pheno.F, pheno.col = "Alb15WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat15wk.F, snps = GM_snps)
# save(qtl.log.Alb15WK.192.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.F.Rdata")
# #Males
# #log Albumin 6WKS
# qtl.log.Alb6WK.192.M <- scanone( pheno = pheno.M, pheno.col = "Alb6WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat6wk.M, snps = GM_snps)
# save(qtl.log.Alb6WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.M.Rdata")
# #log Albumin 10WKS
# qtl.log.Alb10WK.192.M <- scanone( pheno = pheno.M, pheno.col = "Alb10WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat10wk.M, snps = GM_snps)
# save(qtl.log.Alb10WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.M.Rdata")
# #log Albumin 15WKS
# qtl.log.Alb15WK.192.M <- scanone( pheno = pheno.M, pheno.col = "Alb15WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat15wk.M, snps = GM_snps)
# save(qtl.log.Alb15WK.192.M, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.M.Rdata")

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

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
## Additional QTL plotting:
#	Plotting QTL as .pdf without X chromosome
#	Plotting allele effect of baysian interval
#	Plotting average of founder effect by founder
## Yuka Takemon
## 12/05/16

## Setup
library(DOQTL)
library(ggplot2)
library(reshape2)
library(knitr)
setwd("/hpcdata/ytakemon/Col4a5xDO")
#	load files
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15WK.192.Rdata")
#	Create threshold by permutations
thr.1000.qtl.logAlb.6wk <- get.sig.thr( perms.1000.qtl.log.Alb6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.10wk <- get.sig.thr( perms.1000.qtl.log.Alb10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.15wk <- get.sig.thr( perms.1000.qtl.log.Alb15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

#	If you don't know you bayesian interval
#interval = bayesint(qtl, chr = 2)
#knitr::kable(interval)

## Helper functions (http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions)
#	Summarizes data, handling within-subjects variables by removing inter-subject variability.
#	It will still work if there are no within-S variables.
#	Gives count, un-normed mean, normed mean (with same between-group mean),
#	standard deviation, standard error of the mean, and confidence interval.
#	If there are within-subject variables, calculate adjusted values using method from Morey (2008).
#	data: a data frame.
#	measurevar: the name of a column that contains the variable to be summariezed
#	betweenvars: a vector containing names of columns that are between-subjects variables
#	withinvars: a vector containing names of columns that are within-subjects variables
#	idvar: the name of a column that identifies each subject (or matched subjects)
#	na.rm: a boolean that indicates whether to ignore NA's
#	conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
    FUN=is.factor, FUN.VALUE=logical(1))
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  # Apply correction from Morey (2008) to the standard error and confidence interval  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                           FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
#	Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#	data: a data frame.
#	measurevar: the name of a column that contains the variable to be summariezed
#	groupvars: a vector containing names of columns that contain grouping variables
#	na.rm: a boolean that indicates whether to ignore NA's
#	conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
#	Norms the data within specified groups in a data frame; it normalizes each
#	subject (identified by idvar) so that they have the same mean, within each group
#	specified by betweenvars.
#	data: a data frame.
#	idvar: the name of a column that identifies each subject (or matched subjects)
#	measurevar: the name of a column that contains the variable to be summariezed
#	betweenvars: a vector containing names of columns that are between-subjects variables
#	na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    library(plyr)

    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
     .fun = function(xx, col, na.rm) {
        c(subjMean = mean(xx[,col], na.rm=na.rm))
      },
      measurevar,
      na.rm
    )

    # Put the subject means with original data
    data <- merge(data, data.subjMean)

    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                               mean(data[,measurevar], na.rm=na.rm)

    # Remove this subject mean column
    data$subjMean <- NULL

    return(data)
}
## Alb6wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb6WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb6wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.6WK QTL perms.1000.noX")
dev.off()

## Alb10wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb10WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb10wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.10wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.10WK QTL perms.1000.noX")
dev.off()

## Alb10wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb15WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb15wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.15wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.15WK QTL perms.1000.noX")
dev.off()

## Plotting narrower allele effect Alb 10wk
#	subset qtl object to region of interest
# LOD
lodA <- qtl$lod$A
lodA <- lodA[lodA$chr == 2,]
lodA <- lodA[lodA$pos > 101.7, ]
lodA <- lodA[lodA$pos < 113.7, ]
# Coef
coefA <- qtl$coef$A #matrix
coefA <- as.data.frame(coefA)
coefA <- coefA[rownames(lodA),]
# replace
qtl$lod$A <- lodA
qtl$coef$A <- coefA
#	plot allele effect
qtl$coef$A[abs(qtl$coef$A) > 2 ] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr2.bayesian.int.pdf", width = 10.0, height = 7.5)
coefplot(qtl, chr = 2, main = "Chr 2 Alb10WK_log Allele Effect Plot @ bayesian interval")
dev.off()

## Ploting average founder effect of each strain
#	subset out coefficent data
data <- qtl$coef$A
colnames(data)[1] <- "A"
data$sex <- NULL
data$creat10wk <- NULL
#	Center the coefficent values
data[,2:ncol(data)] <- data[,2:ncol(data)] + data[,1]
data <- data - rowMeans(data)
#	reshpe melt dataframe for ggplot
Founder_names <- c("A", "B", "C","D","E","F","G","H")
ggdata <- melt(data, variable.name = "Founders", value.name = "Effect")
#	Get standard error
ggdata_SE <- summarySEwithin(ggdata, measurevar = "Effect", withinvars = "Founders")

ggplot(ggdata, aes(Founders, Effect)) +
		geom_point()+
		geom_errorbar(aes(ymin = Effecct - se, ymax = Effect + se), width = 0.5)


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_founder_effect_QTLAlb6wk_BayInt_Chr2.pdf", width = 10.0, height = 7.5)
ggplot(ggdata_SE, aes(Founders, Effect)) +
	geom_point()+
	geom_errorbar(aes(ymin = Effect - se, ymax = Effect + se), width = 0.5) +
	labs( title = "Average founder effect within bayesian interval of Alb10wk QTL Chr 2", x = "DO Founders", y = "Founder Effect") +
	theme( plot.title = element_text(hjust = 0.5))
dev.off()
