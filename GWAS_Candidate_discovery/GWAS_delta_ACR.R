#Plotting QTL of the delta between ACR at 15 and 6WKs.
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
pheno$delta_ACR15_6 <- pheno$ACR15WK - pheno$ACR6WK
pheno[pheno ==  -Inf] = NA
options(na.action = "na.pass")
pheno <- pheno[,c("MouseID", "Sex", "ACR6WK", "ACR15WK", "delta_ACR15_6")]

#sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"
sex.covar <- as.matrix(sex.covar)

#GWAS
gwas_delta_ACR15_6 <- scanone.assoc(
	pheno = pheno,
	pheno.col = "delta_ACR15_6",
	probs = best.genoprobs.192,
	K = K_GS,
	addcovar = sex.covar,
	markers = GM_snps,
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz",
	ncl = 3) #number of cores to use
save(gwas_delta_ACR15_6, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.delta_ACR15_6.Rdata")

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS_delta_ACR15_6.pdf", width = 11, height = 7)
plot(gwas_delta_ACR15_6, main = "Col4a5xDO delta ACR15_6wk GWAS map")
dev.off()

#####################################################################################################################
##Coefficient plots using QTL Rdata
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl_delta_ACR15_6.192.Rdata")

#Chr1 coef plot
qtl <- qtl.deltaACR15_6.192
qtl$coef$A[abs(qtl$coef$A) > 5000] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr1.pdf", width = 12, height = 6)
coefplot(qtl, chr = 1, main = "Allele Effect for delta_ACR15_6 @ Chr 1")
dev.off()

#Chr2 coef plot
qtl <- qtl.deltaACR15_6.192
qtl$coef$A[abs(qtl$coef$A) > 2000] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr2.pdf", width = 12, height = 6)
coefplot(qtl, chr = 2, main = "Allele Effect for delta_ACR15_6 @ Chr 2")
dev.off()

#Chr5 coef plot
qtl <- qtl.deltaACR15_6.192
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr5.pdf", width = 12, height = 6)
coefplot(qtl, chr = 5, main = "Allele Effect for delta_ACR15_6 @ Chr 5")
dev.off()

#Chr8 coef plot
qtl <- qtl.deltaACR15_6.192
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr8.pdf", width = 12, height = 6)
coefplot(qtl, chr = 8, main = "Allele Effect for delta_ACR15_6 @ Chr 8")
dev.off()

#Chr13 coef plot
qtl <- qtl.deltaACR15_6.192
qtl$coef$A[abs(qtl$coef$A) > 3000] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr13.pdf", width = 12, height = 6)
coefplot(qtl, chr = 13, main = "Allele Effect for delta_ACR15_6 @ Chr 13")
dev.off()

#Chr14 coef plot
qtl <- qtl.deltaACR15_6.192
#qtl$coef$A[abs(qtl$coef$A) > 3000] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr14.pdf", width = 12, height = 6)
coefplot(qtl, chr = 14, main = "Allele Effect for delta_ACR15_6 @ Chr 14")
dev.off()

#Chr17 coef plot
qtl <- qtl.deltaACR15_6.192
qtl$coef$A[abs(qtl$coef$A) > 1500] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr17.pdf", width = 12, height = 6)
coefplot(qtl, chr = 17, main = "Allele Effect for delta_ACR15_6 @ Chr 17")
dev.off()

#Chr19 coef plot
qtl <- qtl.deltaACR15_6.192
#qtl$coef$A[abs(qtl$coef$A) > 1500] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef_delta_ACR15_6_chr19.pdf", width = 12, height = 6)
coefplot(qtl, chr = 19, main = "Allele Effect for delta_ACR15_6 @ Chr 19")
dev.off()

#####################################################################################################################
# Candidate genes
#Chr1 candidate gene proximal
chr <- 1
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 35,
										end = 45,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr1_proximal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr1 candidate gene distal
chr <- 1
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 65,
										end = 75,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr1_distal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr2 candidate gene proximal
chr <- 2
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 30,
										end = 45,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr2_proximal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3.8, show.sdps = TRUE)
dev.off()

#Chr2 candidate gene distal
chr <- 2
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 78,
										end = 92,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr2_distal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 5, show.sdps = TRUE)
dev.off()

#Chr5 candidate gene proximal
chr <- 5
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 45,
										end = 55,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr5_proximal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr5 candidate gene distal
chr <- 5
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 64,
										end = 75,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr5_distal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 5, show.sdps = TRUE)
dev.off()

#Chr8 candidate gene
chr <- 8
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 66,
										end = 72,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr8.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3.5, show.sdps = TRUE)
dev.off()

#Chr13 candidate gene
chr <- 13
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 38,
										end = 42,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr13.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3.5, show.sdps = TRUE)
dev.off()

#Chr14 candidate gene proximal
chr <- 14
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 50,
										end = 60,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr14_proximal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr14 candidate gene medial
chr <- 14
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 65,
										end = 75,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr14_medial.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr14 candidate gene distal
chr <- 14
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 75,
										end = 90,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr14_distal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr17 candidate gene poximal
chr <- 17
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 24,
										end = 30,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr17_poximal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 4, show.sdps = TRUE)
dev.off()

#Chr17 candidate gene distal
chr <- 17
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 65,
										end = 75,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr17_distal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 4, show.sdps = TRUE)
dev.off()

#Chr19 candidate gene proximal
chr <- 19
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 20,
										end = 30,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr19_proximal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

#Chr19 candidate gene distal
chr <- 19
genes <- assoc.map(	pheno = pheno,
										pheno.col = "delta_ACR15_6",
										probs = best.genoprobs.192,
										K = K_GS[[chr]],
										addcovar = sex.covar,
										snps = GM_snps,
										chr = chr,
										start = 45,
										end = 55,
										output = "p-value")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Candidate_delta_ACR15_6_chr19_distal.pdf", width = 12, height = 6)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()
