library(DOQTL)
library(abind)
library(knitr)

dir <- "/projects/marralab/ytakemon_prj/Col4a5/"

#load sample
load(paste0(dir,"Data/consolidated/best.genoprobs.192.Rdata"))
load(paste0(dir,"Data/Col4a5xDO_192data_YT.Rdata"))
load(paste0(dir,"Data/Col4a5xDO_GWAS_YT.Rdata"))
load(paste0(dir,"Data/Col4a5xDO_GWASperm1000_YT.Rdata"))
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#read and clean up phenotype data
pheno <- read.delim(paste0(dir,"Data/consolidated/Phenotype/1415_master_pheno.txt"), sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA # remove negative values after log
pheno[pheno ==  -Inf] = NA # remove non- numbers
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

# Just sex as covariate
sex.covar <- as.data.frame(Covar[,"sex"])
rownames(sex.covar) <- rownames(Covar)
colnames(sex.covar) <- "Sex"
sex.creat.covar6 <- sex.covar
sex.creat.covar6$Creat6WK <- Pheno$Creat6WK
sex.creat.covar10 <- sex.covar
sex.creat.covar10$Creat10WK <- Pheno$Creat10WK
sex.creat.covar15 <- sex.covar
sex.creat.covar15$Creat15WK <- Pheno$Creat15WK

# Permutation threshold
gwas.Alb6.perm <- get.sig.thr( -log10(GWAS.Alb6WK.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.Alb10.perm <- get.sig.thr( -log10(GWAS.Alb10WK.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.Alb15.perm <- get.sig.thr( -log10(GWAS.Alb15WK.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
gwas.gfr.perm <- get.sig.thr(-log10(GWAS.GFR.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

# GWAS plots
#Alb6
png(paste0(dir,"Results/consolidated/best.genoprobs.192.Rdata")"./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.Alb6wk.png", width = 2000, height = 1000, res = 100)
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

# Alb6 chr2
chr <- 2
gwas <- GWAS.Alb6WK
pdf(paste0(dir,"/Results/GWAS.Alb6WK.Chr",chr,".pdf"), width = 12, height = 5)
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb6WK"), sig.thr = gwas.Alb6.perm["0.05"], sig.col = "red", ylim = c(0 ,max(gwas.Alb10.perm)))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb6WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar6,
	snps = Snps,
	chr = chr,
	start = 112,
	end = 125,
	output = "p-value")
pdf(paste0(dir,"/Results/GWAS.Alb6WK.Chr",chr,".candidates.pdf"), width = 15, height = 7)
assoc.plot(genes, show.sdps = TRUE)
dev.off()

# Alb6 chr11
chr <- 11
gwas <- GWAS.Alb6WK
pdf(paste0(dir,"/Results/GWAS.Alb6WK.Chr",chr,".pdf"), width = 12, height = 5)
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb6WK"), sig.thr = gwas.Alb6.perm["0.05"], sig.col = "red", ylim = c(0 ,max(gwas.Alb10.perm)))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb6WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar6,
	snps = Snps,
	chr = chr,
	start = 85,
	end = 92,
	output = "p-value")
pdf(paste0(dir,"/Results/GWAS.Alb6WK.Chr",chr,".candidates.pdf"), width = 15, height = 7)
assoc.plot(genes, thr = 6, show.sdps = TRUE)
dev.off()

# Alb10 chr2
chr <- 2
gwas <- GWAS.Alb10WK
pdf(paste0(dir,"/Results/GWAS.Alb10WK.Chr",chr,".pdf"), width = 12, height = 5)
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb10WK"), sig.thr = gwas.Alb10.perm["0.05"], sig.col = "red", ylim = c(0 ,max(gwas.Alb10.perm)))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb10WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar10,
	snps = Snps,
	chr = chr,
	start = 112,
	end = 125,
	output = "p-value")
pdf(paste0(dir,"/Results/GWAS.Alb10WK.Chr",chr,".candidates.pdf"), width = 15, height = 7)
assoc.plot(genes, show.sdps = TRUE)
dev.off()

chr <- 13
pdf(paste0(dir,"/Results/GWAS.Alb10WK.Chr",chr,".pdf"), width = 12, height = 5)
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb10WK"), sig.thr = gwas.Alb10.perm["0.05"], sig.col = "red", ylim = c(0 ,max(gwas.Alb10.perm)))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb10WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar10,
	snps = Snps,
	chr = chr,
	start = 90,
	end = 100,
	output = "p-value")
pdf(paste0(dir,"/Results/GWAS.Alb10WK.Chr",chr,".candidates.pdf"), width = 15, height = 7)
assoc.plot(genes, show.sdps = TRUE)
dev.off()

# Alb15 chr15
chr <- 15
gwas <- GWAS.Alb15WK
pdf(paste0(dir,"/Results/GWAS.Alb15WK.Chr",chr,".pdf"), width = 12, height = 5)
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb15WK"), sig.thr = gwas.Alb15.perm["0.05"], sig.col = "red", ylim = c(0 ,max(gwas.Alb15.perm)))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb15WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar15,
	snps = Snps,
	chr = chr,
	start = 86,
	end = 90,
	output = "p-value")
pdf(paste0(dir,"/Results/GWAS.Alb15WK.Chr",chr,".candidates.pdf"), width = 15, height = 7)
assoc.plot(genes, show.sdps = TRUE)
dev.off()

# GFR
chr <- 19
gwas <- GWAS.log.C2
pdf(paste0(dir,"/Results/GWAS.GFR.Chr",chr,".pdf"), width = 12, height = 5)
plot(gwas , chr = chr, main = paste("Chr", chr, "GFR"), sig.thr = gwas.gfr.perm["0.05"], sig.col = "red", ylim = c(0 ,max(gwas.gfr.perm)))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="C2_log",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.covar,
	snps = Snps,
	chr = chr,
	start = 26,
	end = 30,
	output = "p-value")
pdf(paste0(dir,"/Results/GWAS.GFR.Chr",chr,".candidates.pdf"), width = 15, height = 7)
assoc.plot(genes, show.sdps = TRUE)
dev.off()
