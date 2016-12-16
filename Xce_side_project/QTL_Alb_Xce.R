library(DOQTL)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
load("./Phenotype/Xce_allele_MF.Rdata") #Xce allele of both male and female 
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)


#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(Xce_allele),] #subset pheno to match 192 samples
pheno$Xce_winner <- Xce_allele$Xce_winner

#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA

pheno$C2_log <- log(pheno$C2) 

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

#keep only interested columns
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Creat6WK_log","Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log", "Xce_winner")]

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#add Xce on or off covariate:
Xce_OF <- model.matrix(~0+Xce_winner, data = pheno)
colnames(Xce_OF)[2] <- "Xce_on_off"
Xce_OF <- as.data.frame(Xce_OF)
sex.covar$Xce_on_off <- Xce_OF$Xce_on_off

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

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")

#QTL mapping of Alb with sex as covariate at 3 different time points
#Both sexes
#GFR
qtl.log.GFR.Xce.192 <- scanone( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#log ACR6WKS
qtl.log.Alb6WK.192 <- scanone( pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat6wk, snps = GM_snps)
#log Alb10WKS
qtl.log.Alb10WK.192 <- scanone( pheno = pheno, pheno.col = "Alb10WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat10wk, snps = GM_snps)
#log Alb15WKS
qtl.log.Alb15WK.192 <- scanone( pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat15wk, snps = GM_snps)

#save qtl files:
save(qtl.log.GFR.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.GFR.Xce.192.Rdata")
qtl.log.Alb6WK.Xce.192 <- qtl.log.Alb6WK.192
save(qtl.log.Alb6WK.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/save(qtl.log.Alb6WK.Xce.192.Rdata")
qtl.log.Alb10WK.Xce.192 <- qtl.log.Alb10WK.192
save(qtl.log.Alb10WK.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/save(qtl.log.Alb10WK.Xce.192.Rdata")
qtl.log.Alb15WK.Xce.192 <- qtl.log.Alb15WK.192
save(qtl.log.Alb15WK.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/save(qtl.log.Alb15WK.Xce.192.Rdata")

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.GFR.Xce.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.GFR.Xce.192, main = "Col4a5xDO log GFR QTL w/Xce")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb6wk.Xce.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.Xce.192, main = "Col4a5xDO log.Alb.6WK QTL w/Xce")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb10wk.Xce.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.Xce.192, main = "Col4a5xDO log.Alb.10WK QTL w/Xce")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb15wk.Xce.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.Xce.192, main = "Col4a5xDO log.Alb.15WK QTL w/Xce")
dev.off()

#PERMUTATIONS
#log GFR
perms.1000.qtl.log.GFR.Xce.192  <- scanone.perm( pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.GFR.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.GFR.Xce.192.Rdata" )
#log Alb6
perms.1000.qtl.log.Alb6WK.Xce.192  <- scanone.perm( pheno = pheno, pheno.col = "Alb6WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.Alb6WK.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.Xce.192.Rdata" )
#log Alb10
perms.1000.qtl.log.Alb10WK.Xce.192  <- scanone.perm( pheno = pheno, pheno.col = "Alb10WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.Alb10WK.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.Xce.192.Rdata" )
#log Alb15
perms.1000.qtl.log.Alb15WK.Xce.192  <- scanone.perm( pheno = pheno, pheno.col = "Alb15WK_log", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
save(perms.1000.qtl.log.Alb15WK.Xce.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15.Xce.192.Rdata" )



#########################start here when perms are done
#load qtls
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.GFR.Xce.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/save(qtl.log.Alb6WK.Xce.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/save(qtl.log.Alb10WK.Xce.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/save(qtl.log.Alb15WK.Xce.192.Rdata")
#load perms
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.GFR.Xce.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.Xce.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.Xce.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15WK.Xce.192.Rdata")

#Assign threshold
#threshold: 1000
thr.1000.qtl.logGFR <- get.sig.thr( perms.1000.qtl.log.GFR.Xce.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)
thr.1000.qtl.logAlb.6wk <- get.sig.thr( perms.1000.qtl.log.Alb6WK.Xce.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)
thr.1000.qtl.logAlb.10wk <- get.sig.thr( perms.1000.qtl.log.Alb10WK.Xce.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)
thr.1000.qtl.logAlb.15wk <- get.sig.thr( perms.1000.qtl.log.Alb15WK.Xce.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)

#plot with perm1000
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.GFR.Xce.perm1000.X_T.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.GFR.Xce.192, sig.thr = thr.1000.qtl.logGFR, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log GFR QTL w/Xce perm1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb6wk.Xce.perm1000.X_T.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.Xce.192, sig.thr = thr.1000.qtl.logAlb.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.6WK QTL w/Xce perm1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb10wk.Xce.perm1000.X_T.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.Xce.192, sig.thr = thr.1000.qtl.logAlb.10wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.10WK QTL w/Xce perm1000")
dev.off() 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb15wk.Xce.perm1000.X_T.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.Xce.192, sig.thr = thr.1000.qtl.logAlb.15wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.15WK QTL w/Xce perm1000")
dev.off()


















