########################################################################################################################
## Compare Fmn1 and Gremline expression
## Author: Yuka Takemon
## Date created: 12/16/17

## OBJECTIVE
## Ron's advisor on the Col4a5 grand pointed out that there is controversy surrounding Fmn1 and Gremlin. Part of gremline's 
## regulatory region is location in the intron of Fmn1, and many previously published experiemtns show that if Fmn1 is 
## manipulated so will Gremlin. Thus leading to uncertainty becase which of the two leads to the phenotype. We can to 
## check the correlation between Fmn1 nad Gremlin in our data set tp make sure we can avoid such uncertainty and confusion 
## in the future. 

## sessionInfo()
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## locale:
## [1] C
## attached base packages:
## [1] grid      parallel  stats4    stats     graphics  grDevices utils
## [8] datasets  methods   base
## other attached packages:
## [1] ggplot2_2.1.0        DOQTL_1.0.0          AnnotationDbi_1.28.2
## [4] GenomeInfoDb_1.2.5   IRanges_2.0.1        S4Vectors_0.4.0
## [7] Biobase_2.26.0       BiocGenerics_0.12.1  RSQLite_1.0.0
## [10] DBI_0.3.1
## loaded via a namespace (and not attached):
## [1] Biostrings_2.34.1      GenomicRanges_1.18.4   MUGAExampleData_1.0.0
## [4] QTLRel_0.2-14          RCurl_1.95-4.7         RUnit_0.4.30
## [7] Rcpp_0.11.3            Rsamtools_1.18.3       XML_3.98-1.3
## [10] XVector_0.6.0          annotate_1.44.0        annotationTools_1.40.0
## [13] biomaRt_2.22.0         bitops_1.0-6           colorspace_1.2-6
## [16] corpcor_1.6.8          gdata_2.17.0           gtable_0.1.2
## [19] gtools_3.5.0           hwriter_1.3.2          mclust_5.1
## [22] munsell_0.4.2          org.Hs.eg.db_3.0.0     org.Mm.eg.db_3.0.0
## [25] plyr_1.8.3             scales_0.3.0           tools_3.1.1
## [28] xtable_1.8-0           zlibbioc_1.12.0

## List of data saved from this script (time savers for reanalysis)
## Rdata:
## Tables:
## Plots:
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_Gremlin1_cor.pdf", width = 10.0, height = 7.5)
########################################################################################################################
library(DOQTL)
library(ggplot2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")

RNA_seqZ <- t(RNA_seq)
RNA_seqZ <- apply(RNA_seqZ, 2, rankZ)
RNA_seqZ <- t(RNA_seqZ)
RNA_seqZ <- as.data.frame(RNA_seqZ)
#Fmn1 		: ENSMUSG00000044042
#Gremlin1 	: ENSMUSG00000074934

fit <- lm(ENSMUSG00000074934 ~ ENSMUSG00000044042, RNA_seqZ)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot <- ggplot(RNA_seqZ, aes( x = ENSMUSG00000044042, y = ENSMUSG00000074934)) +
	geom_smooth(method = lm) +
	geom_point() +
	annotate( "text", x = 1, y = 1, label = eq, fontface = "bold", size = 3) +
	labs( title = "Fmn1 vs Gremlin1 TPM ", x = "Fmn1 ranked TPM", y = "Gremlin1 ranked TPM") +
	theme( plot.title = element_text(hjust = 0.5))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_Gremlin1_cor.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()


