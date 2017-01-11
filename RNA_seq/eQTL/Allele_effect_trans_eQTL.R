################################################################################################################################
## Allele Effect Plot of trans_eQTL genes significant at Fmn1
## Yuka Takemon
## Created: 12/20/16
## Last modified: 12/20/16

## OBJECTIVE
## For each eQTL that had a significant peak at Fmn1 ("JAX00098823" or "rs27466023") we will create an allele effect plot of Chr 2.
## Starting with the 35 eQTL that passed the 0.05 significance threshold. 
## The allele peak will be for the greater Chr 2 region, however a close up can be made later on
## if necessary.

## sessionInfo()
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## locale:
## [1] C
## attached base packages:
## [1] grid      parallel  stats4    stats     graphics  grDevices utils
## [8] datasets  methods   base
## other attached packages:
## [1] knitr_1.11           DOQTL_1.0.0          AnnotationDbi_1.28.2
## [4] GenomeInfoDb_1.2.5   IRanges_2.0.1        S4Vectors_0.4.0
## [7] Biobase_2.26.0       BiocGenerics_0.12.1  RSQLite_1.0.0
## [10] DBI_0.3.1
## loaded via a namespace (and not attached):
## [1] Biostrings_2.34.1      GenomicRanges_1.18.4   MUGAExampleData_1.0.0
## [4] QTLRel_0.2-14          RCurl_1.95-4.7         RUnit_0.4.30
## [7] Rcpp_0.11.3            Rsamtools_1.18.3       XML_3.98-1.3
## [10] XVector_0.6.0          annotate_1.44.0        annotationTools_1.40.0
## [13] biomaRt_2.22.0         bitops_1.0-6           colorspace_1.2-6
## [16] corpcor_1.6.8          digest_0.6.8           gdata_2.17.0
## [19] gtable_0.1.2           gtools_3.5.0           hwriter_1.3.2
## [22] mclust_5.1             munsell_0.4.2          org.Hs.eg.db_3.0.0
## [25] org.Mm.eg.db_3.0.0     plyr_1.8.3             scales_0.3.0
## [28] tools_3.1.1            xtable_1.8-0           zlibbioc_1.12.0

## List of data saved from this script (time savers for reanalysis)
## Rdata:
## Tables:
## Plots:
## pdf( paste(plot_dir, "ENSMUSG00000023236.Scg5.chr2.candidates.pdf"), width = 10.0, height = 7.5)
## pdf( paste(plot_dir, "ENSMUSG00000041219.Arhgap11a.chr2.candidates.pdf"), width = 10.0, height = 7.5)
################################################################################################################################

#	load library
library(DOQTL)
library(knitr)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#	Input data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Fmn1_trans_eQTL.Rdata")
Fmn1_trans_eQTL_Sig.0.05 <- Fmn1_trans_eQTL[Fmn1_trans_eQTL$LOD_score_at_JAX00098823 > 8.104061,]

#	Quick directories
qtl_dir <- "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/"
plot_dir <- "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/"

#	Name reference
Ensemble_ID <- rownames(Fmn1_trans_eQTL_Sig.0.05)
MGI_ID <- Fmn1_trans_eQTL_Sig.0.05$MGI_ID

#	Fist run through
for (i in 1:length(Ensemble_ID)){
	sample_file <- paste(qtl_dir, Ensemble_ID[i], ".", MGI_ID[i], ".eQTL.rankZ.tpm.Rdata", sep = "")
	sample_obj <- paste(Ensemble_ID[i],".",MGI_ID[i],".eQTL", sep = "")
	load(sample_file)
	temp <- get(ls(pattern = sample_obj)) 
	pdf(paste(plot_dir, Ensemble_ID[i], ".", MGI_ID[i],".eQTL.allele.effect.pdf", sep = ""), width = 10.0, heigh = 7.5)
	coefplot(temp, chr = 2, main = paste("Chr2","allele effect plot of", Ensemble_ID[i], MGI_ID[i], "eQTL", sep = " "))
	dev.off()
}
#	Ugh! so much noise! 

## Already able to rule out noisy ones that did not have interval in region of interest
#	subset genes with interval within region of interest, manually
Interest <- Fmn1_trans_eQTL_Sig.0.05[c(3,5,11,12,13,14,17,20,21,23,28,29,30),]
Ensemble_ID <- rownames(Interest)
MGI_ID <- Interest$MGI_ID

#	second run through with subset.
for (i in 1:length(Ensemble_ID)){
	sample_file <- paste(qtl_dir, Ensemble_ID[i], ".", MGI_ID[i], ".eQTL.rankZ.tpm.Rdata", sep = "")
	sample_obj <- paste(Ensemble_ID[i],".",MGI_ID[i],".eQTL", sep = "")
	load(sample_file)
	temp <- get(ls(pattern = sample_obj)) 
	temp$coef$A[abs(temp$coef$A) > 0.5 ] = 0 
	pdf(paste(plot_dir, Ensemble_ID[i], ".", MGI_ID[i],".eQTL.allele.effect.pdf", sep = ""), width = 10.0, heigh = 7.5)
	coefplot(temp, chr = 2, main = paste("Chr2","allele effect plot of", Ensemble_ID[i], MGI_ID[i], "eQTL", sep = " "))
	dev.off()
}


#	look up bayesian interval for each in subset and make sure JAX00098823 falls within interval (113.4831)
for (i in 1:length(Ensemble_ID)){
	sample_obj <- paste(Ensemble_ID[i],".",MGI_ID[i],".eQTL", sep = "")
	print(sample_obj)
	sample_file <- paste(qtl_dir, Ensemble_ID[i], ".", MGI_ID[i], ".eQTL.rankZ.tpm.Rdata", sep = "")
	load(sample_file)
	temp <- get(ls(pattern = sample_obj)) 
	interval <- bayesint(temp , chr = 2)
	print(knitr::kable(interval))
}


#View region in detail for Scg5 and Arhgap11a
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

interval <- bayesint(ENSMUSG00000023236.Scg5.eQTL, chr  = 2)
kable(interval)
chr <- 2

Scg5_chr2 <- assoc.map(
	pheno = RNA_seqZ,
	pheno.col = "ENSMUSG00000023236",
	probs = best.genoprobs.192,
	K = K_GS[[chr]],
	addcovar = sex.covar,
	snps = GM_snps, 
	chr = chr, 
	start = 113,
	end = 114,
	output = "p-value")

Arh_chr2 <- assoc.map(
	pheno = RNA_seqZ,
	pheno.col = "ENSMUSG00000041219",
	probs = best.genoprobs.192,
	K = K_GS[[chr]],
	addcovar = sex.covar,
	snps = GM_snps, 
	chr = chr, 
	start = 113,
	end = 114,
	output = "p-value")

pdf( paste(plot_dir, "ENSMUSG00000023236",".Scg5",".chr2.candidates.pdf"), width = 10.0, height = 7.5)
assoc.plot(Scg5_chr2, thr = 10, show.sdps = TRUE)
dev.off()

pdf( paste(plot_dir, "ENSMUSG00000041219",".Arhgap11a",".chr2.candidates.pdf"), width = 10.0, height = 7.5)
assoc.plot(Arh_chr2, thr = 10, show.sdps = TRUE)
dev.off()

