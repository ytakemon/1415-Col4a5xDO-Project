# Yuka Takemon
# 08/29/18
# Create QTL and GWAS maps for publication
library(DOQTL)
library(biomaRt)
#setwd("/hpcdata/ytakemon/Col4a5xDO")
# set up biomaRt
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

#	load files
load("./Col4a5xDO_192data_YT.Rdata")
load("./Col4a5xDO_QTL_YT.Rdata")
load("./Col4a5xDO_QTLperm1000_YT.Rdata")
load("./Col4a5xDO_GWAS_YT.Rdata")
load("./Col4a5xDO_GWASperm1000_YT.Rdata")

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

# GFR--------------------------------------------------------------------------
## Threshold : QTL
# Calcualte QTL threshold
thr <- get.sig.thr( perms.1000.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
# show threshold
thr
#    0.05      0.1     0.63
#7.340369 6.945842 4.659221

#	remove X chr from qtl object
qtl <- qtl.GFR.log.C2.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL

# Plot QTL: GFR
pdf("../Results/QTL.log.C2.GFR.noX.pdf", width = 12, height = 5)
plot(qtl, sig.thr = thr,
  sig.col = c("red", "orange", "chartreuse"),
  main = "Col4a5xDO log.C2.GFR QTL perms.1000.noX")
dev.off()

# Generate gene list in invterval:
bayes <- bayesint(qtl, 19, 0.95)
#> bayes
#                 marker chr      pos       cM perc.var      lrs      lod
#UNC30136851 UNC30136851  19 26.98658 14.99423 20.10076 27.82608 6.042356
#JAX00474835 JAX00474835  19 27.61228 15.38569 24.29844 34.51806 7.495501
#UNCHS047487 UNCHS047487  19 28.04451 15.68458 16.49212 22.34842 4.852898
#                       p neg.log10.p
#UNC30136851 2.364325e-04    3.626293
#JAX00474835 1.377513e-05    4.860904
#UNCHS047487 2.211717e-03    2.655270
chr <- 19
interval <- bayes$pos * 1e6
bayes_geneList <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"),
  filters = c("chromosome_name", "start", "end"),
  values = list(chr, interval[1], interval[3]),
  mart = ensembl)
write.csv(bayes_geneList, "../Results/GFR_chr19_bayesint_geneList.csv",
  row.names = FALSE,
  quote = FALSE)

# Plot GWAS : GFR
thr <- get.sig.thr(-log10(GWAS.GFR.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
# plot GWAS: GFR
pdf("../Results/GWAS.log.C2.GFR.noX.pdf", width = 12, height = 5) # all GWAS
plot(GWAS.log.C2,
  ylim = c(0 ,max(thr)),
  sig.thr = thr,
  sig.col = c("red", "orange", "chartreuse"),
  main = "Col4a5xDO GWAS of log C2 GFR")
dev.off()
pdf("../Results/GWAS.log.C2.GFR.chr19.pdf", width = 12, height = 5) # chr 19 GWAS
plot(GWAS.log.C2 , chr = 19, main = "Chr 19 log.C2.GFR")
dev.off()
pdf("../Results/Coef.log.C2.GFR.chr19.pdf", width = 12, height = 5) # Chr 19 coef plot
coefplot(qtl , chr = 19, main = "Chr 19 log.C2.GFR")
dev.off()

# GFR candiate gene
chr19.genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="C2_log",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.covar,
	snps = Snps,
	chr = 19,
	start = bayes[1,3],
	end = bayes[3,3],
	output = "p-value")
pdf("../Results/GWAS.log.C2.GFR.Chr19.candiates.pdf", width = 15, height = 7)
assoc.plot(chr19.genes, thr = 5, show.sdps = TRUE)
dev.off()

# Alb 6WK-----------------------------------------------------------------------

# Threshold
thr <- get.sig.thr( perms.1000.qtl.log.Alb6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#> thr
#    0.05      0.1     0.63
#7.305748 6.890369 4.756593

#	remove X chr from qtl object
qtl <- qtl.log.Alb6WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL

# Plot QTL: Alb6
pdf("../Results/QTL.log.Alb6WK.noX.pdf", width = 12, height = 5)
plot(qtl, sig.thr = thr,
    sig.col = c("red", "orange", "chartreuse"),
    main = "Col4a5xDO log.Alb6Wk QTL perms.1000.noX")
dev.off()

# Generate gene list in invterval:
chr <- 2
bayes <- bayesint(qtl, chr, 0.95)
#> bayes
#                 marker chr      pos       cM  perc.var      lrs      lod
#JAX00499427 JAX00499427   2 109.1249 49.54448 10.943882 19.23998 4.177908
#JAX00098823 JAX00098823   2 113.4831 50.68170 16.654300 30.24075 6.566695
#UNC4549601   UNC4549601   2 174.0296 87.67960  9.869729 17.24974 3.745734
#                       p neg.log10.p
#JAX00499427 7.467908e-03    2.126801
#JAX00098823 8.577316e-05    4.066649
#UNC4549601  1.585487e-02    1.799837

interval <- bayes$pos * 1e6
bayes_geneList <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"),
  filters = c("chromosome_name", "start", "end"),
  values = list(chr, interval[1], interval[3]),
  mart = ensembl)
write.csv(bayes_geneList, "../Results/Alb6WK_chr2_bayesint_geneList.csv", row.names = FALSE, quote = FALSE)

# Plot GWAS : Alb6WK
thr <- get.sig.thr(-log10(GWAS.Alb6WK.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
# plot GWAS: Alb6WK
gwas <- GWAS.Alb6WK
pdf("../Results/GWAS.log.Alb6WK.noX.pdf", width = 12, height = 5) # all GWAS
plot(gwas,
  ylim = c(0 ,max(thr)),
  sig.thr = thr,
  sig.col = c("red", "orange", "chartreuse"),
  main = "Col4a5xDO GWAS of log Alb6WK")
dev.off()
pdf(paste0("../Results/GWAS.log.Alb6WK.chr",chr,".pdf"), width = 12, height = 5) # chr 19 GWAS
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb6WK"))
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 2 ] = 0
pdf(paste0("../Results/Coef.log.Alb6WK.chr",chr,".pdf"), width = 12, height = 5) # Chr 19 coef plot
coefplot(qtl , chr = chr, main = paste("Chr", chr, "log.Alb6WK"))
dev.off()

# Alb6WK candiate gene
genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb6WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar6,
	snps = Snps,
	chr = chr,
	start = 110, # this is within the bayes interval just narrower
	end = 114,
	output = "p-value")
pdf(paste0("../Results/GWAS.log.Alb6WK.Chr",chr,".candiates.pdf"), width = 15, height = 7)
assoc.plot(genes, thr = 3, show.sdps = TRUE)
dev.off()

# Alb 10WK-----------------------------------------------------------------------

# Threshold
thr <- get.sig.thr( perms.1000.qtl.log.Alb10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#> thr
#    0.05      0.1     0.63
#8.515490 7.776707 5.728013

#	remove X chr from qtl object
qtl <- qtl.log.Alb10WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL

# Plot QTL: Alb10
pdf("../Results/QTL.log.Alb10WK.noX.pdf", width = 12, height = 5)
plot(qtl, sig.thr = thr,
    sig.col = c("red", "orange", "chartreuse"),
    main = "Col4a5xDO log.Alb10Wk QTL perms.1000.noX")
dev.off()

# Generate gene list in invterval:
chr <- 2
bayes <- bayesint(qtl, chr, 0.95)
#> bayes
#                 marker chr      pos       cM perc.var      lrs      lod
#UNCHS005963 UNCHS005963   2 101.7360 45.94599 11.89948 20.39737 4.429233
#JAX00098817 JAX00098817   2 113.4217 50.67787 16.48704 29.00710 6.298813
#UNC3776734   UNC3776734   2 113.6730 50.81143 11.95717 20.50284 4.452134
#                       p neg.log10.p
#UNCHS005963 0.0047723814    2.321265
#JAX00098817 0.0001442559    3.840866
#UNC3776734  0.0045800713    2.339128

interval <- bayes$pos * 1e6
bayes_geneList <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"),
  filters = c("chromosome_name", "start", "end"),
  values = list(chr, interval[1], interval[3]),
  mart = ensembl)
write.csv(bayes_geneList, paste0("../Results/Alb10WK_chr",chr,"_bayesint_geneList.csv"), row.names = FALSE, quote = FALSE)

# Plot GWAS : Alb10WK
thr <- get.sig.thr(-log10(GWAS.Alb10WK.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
# plot GWAS: Alb10WK
gwas <- GWAS.Alb10WK
pdf("../Results/GWAS.log.Alb6WK.noX.pdf", width = 12, height = 5) # all GWAS
plot(gwas,
  ylim = c(0 ,max(thr)),
  sig.thr = thr,
  sig.col = c("red", "orange", "chartreuse"),
  main = "Col4a5xDO GWAS of log Alb10WK")
dev.off()
pdf(paste0("../Results/GWAS.log.Alb10WK.chr",chr,".pdf"), width = 12, height = 5) # chr 19 GWAS
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb10WK"))
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 2 ] = 0
pdf(paste0("../Results/Coef.log.Alb10WK.chr",chr,".pdf"), width = 12, height = 5) # Chr 19 coef plot
coefplot(qtl , chr = chr, main = paste("Chr", chr, "log.Alb10WK"))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb10WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar10,
	snps = Snps,
	chr = chr,
	start = 110, # this is within the bayes interval just narrower
	end = 114,
	output = "p-value")
pdf(paste0("../Results/GWAS.log.Alb10WK.Chr",chr,".candiates.pdf"), width = 15, height = 7)
assoc.plot(genes, thr = 4, show.sdps = TRUE)
dev.off()

# Alb 15WK-----------------------------------------------------------------------

# Threshold
thr <- get.sig.thr( perms.1000.qtl.log.Alb15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#> thr
#    0.05      0.1     0.63
#8.515490 7.776707 5.728013

#	remove X chr from qtl object
qtl <- qtl.log.Alb15WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL

# Plot QTL: Alb10
pdf("../Results/QTL.log.Alb15WK.noX.pdf", width = 12, height = 5)
plot(qtl, sig.thr = thr,
    sig.col = c("red", "orange", "chartreuse"),
    main = "Col4a5xDO log.Alb15Wk QTL perms.1000.noX")
dev.off()

# Generate gene list in invterval:
chr <- 15
bayes <- bayesint(qtl, chr, 0.95)
#> bayes
#                 marker chr      pos       cM perc.var      lrs       lod
#UNC26048887 UNC26048887  15 87.44608 37.22123 22.18870 42.39931  9.206894
#UNCHS041212 UNCHS041212  15 87.67108 37.35712 24.60839 47.73815 10.366207
#UNC26057169 UNC26057169  15 88.00011 37.75250 22.22703 42.48259  9.224977
#                       p neg.log10.p
#UNC26048887 4.355811e-07    6.360931
#UNCHS041212 4.005879e-08    7.397302
#UNC26057169 4.197724e-07    6.376986

interval <- bayes$pos * 1e6
bayes_geneList <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"),
  filters = c("chromosome_name", "start", "end"),
  values = list(chr, interval[1], interval[3]),
  mart = ensembl)
write.csv(bayes_geneList, paste0("../Results/Alb15WK_chr",chr,"_bayesint_geneList.csv"), row.names = FALSE, quote = FALSE)

# Plot GWAS : Alb15WK
thr <- get.sig.thr(-log10(GWAS.Alb15WK.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
# plot GWAS: Alb15WK
gwas <- GWAS.Alb15WK
pdf("../Results/GWAS.log.Alb6WK.noX.pdf", width = 12, height = 5) # all GWAS
plot(gwas,
  ylim = c(0 ,max(thr)),
  sig.thr = thr,
  sig.col = c("red", "orange", "chartreuse"),
  main = "Col4a5xDO GWAS of log Alb15WK")
dev.off()
pdf(paste0("../Results/GWAS.log.Alb15WK.chr",chr,".pdf"), width = 12, height = 5) # chr 19 GWAS
plot(gwas , chr = chr, main = paste("Chr", chr, "log.Alb15WK"))
dev.off()
pdf(paste0("../Results/Coef.log.Alb15WK.chr",chr,".pdf"), width = 12, height = 5) # Chr 19 coef plot
coefplot(qtl , chr = chr, main = paste("Chr", chr, "log.Alb15WK"))
dev.off()

genes <- assoc.map(
	pheno = Pheno,
	pheno.col ="Alb15WK",
	probs = Genoprobs,
	K = K[[chr]],
	addcovar = sex.creat.covar15,
	snps = Snps,
	chr = chr,
	start = bayes[1,3],
	end = bayes[3,3],
	output = "p-value")
pdf(paste0("../Results/GWAS.log.Alb15WK.Chr",chr,".candiates.pdf"), width = 15, height = 7)
assoc.plot(genes, thr = 2, show.sdps = TRUE)
dev.off()
