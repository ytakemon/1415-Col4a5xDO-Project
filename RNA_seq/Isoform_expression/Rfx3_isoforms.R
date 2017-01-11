################################################################################################################################
## Col4a5xDO Rfx isoforms
## Yuka Takemon
## Created: 12/09/16
## Last modified: 12/09/16

## OBJECTIVE
## We want to identify which transcript of Rfx3 is expressed in the Col4a5 B6xDO animals.
##
## Gene Ensembl ID: ENSMUSG00000040929
##
## Transcripts according to Ensembl:
## Rfx3_001: ENSMUST00000174850
## Rfx3_002: ENSMUST00000173161
## Rfx3_003: ENSMUST00000046898 
## Rfx3_004: ENSMUST00000172907
## Rfx3_005: ENSMUST00000173863
## Rfx3_006: ENSMUST00000172498
## Rfx3_007: ENSMUST00000174420
## Rfx3_201: ENSMUST00000165566

## sessionInfo()
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## locale:
## [1] C
## attached base packages:
## [1] grid      parallel  stats4    stats     graphics  grDevices utils
## [8] datasets  methods   base
## other attached packages:
## [1] reshape2_1.4.1       knitr_1.11           DOQTL_1.0.0
## [4] AnnotationDbi_1.28.2 GenomeInfoDb_1.2.5   IRanges_2.0.1
## [7] S4Vectors_0.4.0      Biobase_2.26.0       BiocGenerics_0.12.1
## [10] RSQLite_1.0.0        DBI_0.3.1
## loaded via a namespace (and not attached):
## [1] Biostrings_2.34.1      GenomicRanges_1.18.4   MUGAExampleData_1.0.0
## [4] QTLRel_0.2-14          RCurl_1.95-4.7         RUnit_0.4.30
## [7] Rcpp_0.11.3            Rsamtools_1.18.3       XML_3.98-1.3
## [10] XVector_0.6.0          annotate_1.44.0        annotationTools_1.40.0
## [13] biomaRt_2.22.0         bitops_1.0-6           colorspace_1.2-6
## [16] corpcor_1.6.8          digest_0.6.8           gdata_2.17.0
## [19] gtable_0.1.2           gtools_3.5.0           hwriter_1.3.2
## [22] magrittr_1.5           mclust_5.1             munsell_0.4.2
## [25] org.Hs.eg.db_3.0.0     org.Mm.eg.db_3.0.0     plyr_1.8.3
## [28] scales_0.3.0           stringi_1.0-1          stringr_1.0.0
## [31] tools_3.1.1            xtable_1.8-0           zlibbioc_1.12.0

## List of data saved from this script (time savers for reanalysis)
## Rdata:
## Tables:
## Plots:
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3_transcript_tpm.pdf", width = 10.0, height = 7.5)
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)	
################################################################################################################################

#Load libraries
library(DOQTL)
library(ggplot2)
library(reshape2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")
#load data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")
Rfx3_ID <- c("ENSMUST00000174850", "ENSMUST00000173161", "ENSMUST00000046898", "ENSMUST00000172907", "ENSMUST00000173863", "ENSMUST00000172498", 
			"ENSMUST00000174420", "ENSMUST00000165566")
Rfx3_names <- c("Rfx3_001", "Rfx3_002", "Rfx3_003", "Rfx3_004", "Rfx3_005", "Rfx3_006", "Rfx3_007", "Rfx3_201")
Gene_allele <- as.data.frame(Gene_allele)

All_transcript_tpm <- as.data.frame(All_transcript_tpm)
All_transcript_tpm_rankZ <- as.data.frame(All_transcript_tpm_rankZ)

Rfx3_ID_tpm <- All_transcript_tpm[, Rfx3_ID]
Rfx3_ID_ranked_tpm <- All_transcript_tpm_rankZ[, Rfx3_ID]

gg_data <- Rfx3_ID_tpm
colnames(gg_data) <- Rfx3_names
gg_data$allele <- Gene_allele$ENSMUSG00000040929
gg_data <- melt(gg_data)
names(gg_data) <- c("Allele", "Rfx3_transcript_ID", "Values")

ggplot <- ggplot(gg_data, aes( x = Rfx3_transcript_ID, y = Values, fill = Rfx3_transcript_ID)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.05) +
	scale_x_discrete( "Rfx3 transcripts") +
	scale_y_continuous( "transcript tpm") +
	labs( title = "Comparison of Rfx3 transcirpts tpm") +
	theme( plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3_transcript_tpm.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

ggplot <- ggplot(gg_data, aes( x = Rfx3_transcript_ID, y = Values, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.05, position = position_dodge(1)) +
	scale_x_discrete( "Rfx3 transcripts") +
	scale_y_continuous( "transcript tpm") +
	labs( title = "Comparison of Rfx3 transcirpts tpm") +
	theme( plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()
