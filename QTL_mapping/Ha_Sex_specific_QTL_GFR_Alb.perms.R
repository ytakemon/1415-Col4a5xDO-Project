#Individual female and male QTL permutations for GFR and Alb

library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")
args <- commandArgs(trailingOnly = TRUE)

#load files
#genoprobs
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
#GM_grid
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
#kinship
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.F.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.M.Rdata")
#phenotype
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean up phenotype data and calculate logs
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
#Log transform GFR
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2)
#Log transform Alb
pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Alb10WK_log <- log(pheno$Alb10WK)
pheno$Alb15WK_log <- log(pheno$Alb15WK)
#Log transform Creat
pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno$Creat10WK_log <- log(pheno$Creat10WK)
pheno$Creat15WK_log <- log(pheno$Creat15WK)
#final touch up
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#create a subset of phenotype, genoprobs, and covariate by separating female and male samples
best.genoprobs.192 <- best.genoprobs.192[ order( rownames( best.genoprobs.192)),,]
pheno <- pheno[ order( rownames(pheno)),]
#subset male and female phenotype
pheno.F <- subset(pheno, Sex == "F")
pheno.M <- subset(pheno, Sex == "M")
#subset male and female genoprobs object
genoprob.F <- best.genoprobs.192[rownames(pheno.F),,]
genoprob.M <- best.genoprobs.192[rownames(pheno.M),,]
#subset male and female covariates
sex.covar$MouseID <- rownames(sex.covar)
sex.covar <- sex.covar[ order( sex.covar$MouseID),]
sex.covar.F <- sex.covar[rownames(pheno.F),]
sex.covar.M <- sex.covar[rownames(pheno.M),]
#remove extra column
sex.covar$MouseID <- NULL
sex.covar.F$MouseID <- NULL
sex.covar.M$MouseID <- NULL

#create covar with sex and creatinine for Albumin QTL
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

#Create randomized samples for permutation
#index females and males
females <- which(sex.covar.F[,1] == 0)
males <- which(sex.covar.M[,1] == 1)
#randomized samples 
random.F <- sample( rownames(pheno.F)[females])
random.M <- sample( rownames(pheno.M)[males])
#randomize pheno data
rownames(pheno.F) <- random.F
rownames(pheno.M) <- random.M
#randomize sex.covar
rownames(sex.covar.F) <- random.F
rownames(sex.covar.M) <- random.M

#set as matrix to feed into qlt
sex.covar.F <- as.matrix(sex.covar.F)
covar.sex.creat6wk.F <- as.matrix(covar.sex.creat6wk.F)
covar.sex.creat10wk.F <- as.matrix(covar.sex.creat10wk.F)
covar.sex.creat15wk.F <- as.matrix(covar.sex.creat15wk.F)

sex.covar.M <- as.matrix(sex.covar.M)
covar.sex.creat6wk.M <- as.matrix(covar.sex.creat6wk.M)
covar.sex.creat10wk.M <- as.matrix(covar.sex.creat10wk.M)
covar.sex.creat15wk.M <- as.matrix(covar.sex.creat15wk.M)

#run QTLs
#females
qtl_G.F <- scanone(pheno = pheno.F, pheno.col = "C2_log", probs = genoprob.F, K = K_GS.F, addcovar = sex.covar.F, snps = GM_snps)
qtl_A6.F <- scanone(pheno = pheno.F, pheno.col = "Alb6WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat6wk.F, snps = GM_snps)
qtl_A10.F <- scanone(pheno = pheno.F, pheno.col = "Alb10WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat10wk.F, snps = GM_snps)
qtl_A15.F <- scanone(pheno = pheno.F, pheno.col = "Alb15WK_log", probs = genoprob.F, K = K_GS.F, addcovar = covar.sex.creat15wk.F, snps = GM_snps)
#males
qtl_G.M <- scanone(pheno = pheno.M, pheno.col = "C2_log", probs = genoprob.M, K = K_GS.M, addcovar = sex.covar.M, snps = GM_snps)
qtl_A6.M <- scanone(pheno = pheno.M, pheno.col = "Alb6WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat6wk.M, snps = GM_snps)
qtl_A10.M <- scanone(pheno = pheno.M, pheno.col = "Alb10WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat10wk.M, snps = GM_snps)
qtl_A15.M <- scanone(pheno = pheno.M, pheno.col = "Alb15WK_log", probs = genoprob.M, K = K_GS.M, addcovar = covar.sex.creat15wk.M, snps = GM_snps)

#take max LOD
#females
max_lod_G.F <- max(qtl_G.F$lod$A$lod)
max_lod_A6.F <- max(qtl_A6.F$lod$A$lod)
max_lod_A10.F <- max(qtl_A10.F$lod$A$lod)
max_lod_A15.F <- max(qtl_A15.F$lod$A$lod)
#males
max_lod_G.M <- max(qtl_G.M$lod$A$lod)
max_lod_A6.M <- max(qtl_A6.M$lod$A$lod)
max_lod_A10.M <- max(qtl_A10.M$lod$A$lod)
max_lod_A15.M <- max(qtl_A15.M$lod$A$lod)

#write lod result
#female
write(max_lod_G.F, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.F.perm/lod.perm.",args[1],".txt"))
write(max_lod_A6.F, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.F.perm/lod.perm.",args[1],".txt"))
write(max_lod_A10.F, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.F.perm/lod.perm.",args[1],".txt"))
write(max_lod_A15.F, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.F.perm/lod.perm.",args[1],".txt"))
#male
write(max_lod_G.M, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.M.perm/lod.perm.",args[1],".txt"))
write(max_lod_A6.M, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.M.perm/lod.perm.",args[1],".txt"))
write(max_lod_A10.M, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.M.perm/lod.perm.",args[1],".txt"))
write(max_lod_A15.M, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.M.perm/lod.perm.",args[1],".txt"))

#End of loop


############################
#This portion is in linux Cadillac
grep -h . *.txt > 1000.perms.txt
#check number of lines/ permutation in file
wc -l 1000.perms.txt


###########################
#Switch to interactive R processing for plot creating
#load perms and get threshold
perm.G.F <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.F.perm/1000.perms.txt", header = F)
perm.A6.F <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.F.perm/1000.perms.txt", header = F)
perm.A10.F <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.F.perm/1000.perms.txt", header = F)
perm.A15.F <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.F.perm/1000.perms.txt", header = F)

perm.G.M <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.M.perm/1000.perms.txt", header = F)
perm.A6.M <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.M.perm/1000.perms.txt", header = F)
perm.A10.M <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.M.perm/1000.perms.txt", header = F)
perm.A15.M <- read.delim(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.M.perm/1000.perms.txt", header = F)

thr.G.F <- get.sig.thr( perm.G.F , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.A6.F <- get.sig.thr( perm.A6.F , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.A10.F <- get.sig.thr( perm.A10.F , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.A15.F <- get.sig.thr( perm.A15.F , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

thr.G.M <- get.sig.thr( perm.G.M , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.A6.M <- get.sig.thr( perm.A6.M , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.A10.M <- get.sig.thr( perm.A10.M , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.A15.M <- get.sig.thr( perm.A15.M , alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

#load qtl files
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.F.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.F.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.F.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.F.Rdata")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.M.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.M.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.M.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.logsys.Alb15WK.192.M.Rdata")

#plot map
#females
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.GFR.log.C2.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.GFR.log.C2.F, sig.thr = thr.G.F, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Female log GFR C2 QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.log.Alb6.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192.F, sig.thr = thr.A6.F, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Female log Alb6WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.log.Alb10.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192.F, sig.thr = thr.A10.F, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Female log Alb10WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.log.Alb15.F.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192.F, sig.thr = thr.A15.F, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Female log Alb15WK QTL perms.1000")
dev.off()
#males
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.GFR.log.C2.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.GFR.log.C2.M, sig.thr = thr.G.M, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Male log GFR C2 QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.log.Alb6.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb6WK.192.M, sig.thr = thr.A6.M, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Male log Alb6WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.log.Alb10.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb10WK.192.M, sig.thr = thr.A10.M, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Male log Alb10WK QTL perms.1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.perms.1000.log.Alb15.M.png", width = 1500, height = 1000, res = 100)
plot(qtl.log.Alb15WK.192.M, sig.thr = thr.A15.M, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Male log Alb15WK QTL perms.1000")
dev.off()


