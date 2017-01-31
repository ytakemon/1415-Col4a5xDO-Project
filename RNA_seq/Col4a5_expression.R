########################################################################################################################
## Comparing Type IV Collagen expressions
## Author: Yuka Takemon
## Date created: 01/26/16
 
## OBJECTIVE
## We want to plot all TypeIV Collagen expression to make sure that Col4a5 verify KO of Col4a5

##> sessionInfo()
##R version 3.3.1 (2016-06-21)
##Platform: x86_64-pc-linux-gnu (64-bit)
##Running under: CentOS release 6.5 (Final)
##locale:
##[1] C
##attached base packages:
## [1] grid      stats4    parallel  stats     graphics  grDevices utils
## [8] datasets  methods   base
##other attached packages:
## [1] reshape2_1.4.2                     ggplot2_2.2.0
## [3] DOQTL_1.10.0                       VariantAnnotation_1.18.7
## [5] Rsamtools_1.26.1                   SummarizedExperiment_1.2.3
## [7] Biobase_2.34.0                     BSgenome.Mmusculus.UCSC.mm10_1.4.0
## [9] BSgenome_1.40.1                    rtracklayer_1.32.2
##[11] Biostrings_2.42.0                  XVector_0.14.0
##[13] GenomicRanges_1.26.1               GenomeInfoDb_1.10.1
##[15] IRanges_2.8.1                      S4Vectors_0.12.0
##[17] BiocGenerics_0.20.0
##loaded via a namespace (and not attached):
## [1] mclust_5.2              Rcpp_0.12.8             mvtnorm_1.0-5
## [4] lattice_0.20-34         corpcor_1.6.8           class_7.3-14
## [7] gtools_3.5.0            digest_0.6.10           assertthat_0.1
##[10] foreach_1.4.3           plyr_1.8.4              RSQLite_1.0.0
##[13] zlibbioc_1.20.0         GenomicFeatures_1.26.0  lazyeval_0.2.0
##[16] diptest_0.75-7          annotate_1.52.0         gdata_2.17.0
##[19] kernlab_0.9-25          QTLRel_0.2-15           RUnit_0.4.31
##[22] labeling_0.3            BiocParallel_1.8.1      stringr_1.1.0
##[25] RCurl_1.95-4.8          biomaRt_2.30.0          munsell_0.4.3
##[28] nnet_7.3-12             tibble_1.2              codetools_0.2-15
##[31] XML_3.98-1.4            GenomicAlignments_1.8.4 MASS_7.3-45
##[34] regress_1.3-14          bitops_1.0-6            annotationTools_1.48.0
##[37] xtable_1.8-2            gtable_0.2.0            DBI_0.5-1
##[40] magrittr_1.5            scales_0.4.1            stringi_1.1.2
##[43] hwriter_1.3.2           flexmix_2.3-13          doParallel_1.0.10
##[46] robustbase_0.92-6       iterators_1.0.8         tools_3.3.1
##[49] fpc_2.1-10              trimcluster_0.1-2       DEoptimR_1.0-8
##[52] AnnotationDbi_1.36.0    colorspace_1.3-1        rhdf5_2.16.0
##[55] cluster_2.0.5           prabclus_2.2-6          modeltools_0.2-21
########################################################################################################################

library(DOQTL)
library(ggplot2)
library(grid)
library(reshape2)
setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#convert to df format
RNA_seq <- as.data.frame(RNA_seq)
Gene_allele <- as.data.frame(Gene_allele)

#clean pheno data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
#keep only interested columns
pheno <- pheno[,c("MouseID", "Sex")]

# Add all Col4a1-5 tpm and alleles
## ENSEMBL ID
## Col4a1: ENSMUSG00000031502
## Col4a2: ENSMUSG00000031503
## Col4a3: ENSMUSG00000079465
## Col4a4: ENSMUSG00000067158
## Col4a5: ENSMUSG00000031274

#Add A
#Load pheno data
pheno$Col4a1_tpm <- RNA_seq$ENSMUSG00000031502
pheno$Col4a2_tpm <- RNA_seq$ENSMUSG00000031503
pheno$Col4a3_tpm <- RNA_seq$ENSMUSG00000079465
pheno$Col4a4_tpm <- RNA_seq$ENSMUSG00000067158
pheno$Col4a5_tpm <- RNA_seq$ENSMUSG00000031274

#pheno$Col4a1_allele <- Gene_allele$ENSMUSG00000031502
#pheno$Col4a2_allele <- Gene_allele$ENSMUSG00000031503
#pheno$Col4a3_allele <- Gene_allele$ENSMUSG00000079465
#pheno$Col4a4_allele <- Gene_allele$ENSMUSG00000067158
#pheno$Col4a5_allele <- Gene_allele$ENSMUSG00000031274

#Subset Females and males
pheno_F <- pheno[pheno$Sex == "F",]
pheno_M <- pheno[pheno$Sex == "M",]

#plot female, Col3aX expression
ggdata <- melt(pheno_F)
names(ggdata) <- c("MouseID", "Sex", "ColIV", "TPM")
ggplot <- ggplot(ggdata, aes(x = ColIV, y = TPM, fill = ColIV))+
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2.5) +
	scale_x_discrete("Type IV Collagens") +
	scale_y_continuous("Transcript TPM") +
	labs(title = "Comparison of Type IV Collagen TPM expression in Females") +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/TypeIVCollagen_TPM_Female.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

#plot Male, Col3aX expression
ggdata <- melt(pheno_M)
names(ggdata) <- c("MouseID", "Sex", "ColIV", "TPM")
ggplot <- ggplot(ggdata, aes(x = ColIV, y = TPM, fill = ColIV))+
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 4) +
	scale_x_discrete("Type IV Collagens") +
	scale_y_continuous("Transcript TPM") +
	labs(title = "Comparison of Type IV Collagen TPM expression in Males") +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/TypeIVCollagen_TPM_Male.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

#plot Both sex Col3aX expression by sex
ggdata <- melt(pheno)
names(ggdata) <- c("MouseID", "Sex", "ColIV", "TPM")
ggplot <- ggplot(ggdata, aes(x = ColIV, y = TPM, fill = Sex))+
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 4) +
	scale_x_discrete("Type IV Collagens") +
	scale_y_continuous("Transcript TPM") +
	labs(title = "Comparison of Type IV Collagen TPM expression by Sex") +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/TypeIVCollagen_TPM_bySex.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()










