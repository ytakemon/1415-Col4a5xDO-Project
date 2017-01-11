################################################################################################################################
## Col4a5xDO Fmn1 isoform sequence
## Yuka Takemon
## Created: 12/02/16

## OBJECTIVE
## We want to identify which transcript of Fmn1 is expressed in the Col4a5 B6xDO animals.
## Isoforms have not been previously documented in detail, and we are using Chan 1995 to decipher which transcirpts we are seeing 
##
## Probes I-III
## Chan 1995: Probe I-III was constructed by subcloning the 1692 basepair XbaI/HindIII fragment from isoform I (Woychik et al., 1990b) into pBluescript Stratagene).
##
## Probes IV
## Chan 1995: Probe IV is a 1568 basepair EcoRI/NaeI fragment from isoform IV (Jackson- Grusby et al., 1992) in pBluescript.
##
## Probs I-IV
## Chan 1995: Probe I-IV is a 1346 basepair Stu/Spe fragment from isoform I in pBluescript. 
##
## GeneID Ensembl
## Fmn1 ENSMUSG00000044042
##
## TranscriptID Ensembl
## Fmn1-001: ENSMUST00000102547.9 (20 exons)
## Fmn1-002: ENSMUST00000081349.8 (17 exons)
## Fmn1-003: ENSMUST00000110954.6 (4 exons)
## Fmn1-004: ENSMUST00000152255.1 (5 exons)
## Fmn1-005: ENSMUST00000153151.1 (2 exons)
## Fmn1-006: ENSMUST00000145891.7 (6 exons)
## Fmn1-007: ENSMUST00000150510.7 (5 exons)
## Fmn1-008: ENSMUST00000161731.4 (16 exons)
## Fmn1-201: ENSMUST00000099576.8 (19 exons)
## 
## TranscriptID as found in isoform.tpm files
## Fmn1-001: ENSMUST00000102547
## Fmn1-002: ENSMUST00000081349
## Fmn1-003: ENSMUST00000110954
## Fmn1-004: ENSMUST00000152255
## Fmn1-005: ENSMUST00000153151
## Fmn1-006: ENSMUST00000145891
## Fmn1-007: ENSMUST00000150510
## Fmn1-008: ENSMUST00000161731
## Fmn1-201: ENSMUST00000099576

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
## save(All_transcript_tpm, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")
## save(Fmn1_transcript_tpm, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Fmn1_transcript_tpm.Rdata")

## Tables:
## Plots:
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcirpt_tpm_rankZ.png", width = 1500, height = 1000, res =
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcirpt_tpm.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_002.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_005.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_201.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_002.coef.chr2.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_005.coef.chr2.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_201.coef.chr2.png", width = 1500, height = 1000, res = 100)
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcript_rankZ_tpm_by_allele.pdf", width = 10.0, height = 7.5)
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_coding_transcript_tpm.pdf", width = 10.0, height = 7.5)
## pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_coding_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)

################################################################################################################################
library(DOQTL)
library(ggplot2)
library(reshape2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load essential data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
Gene_allele <- as.data.frame(Gene_allele)

#Get directory list for each sample's civet run results
civet_dir <- list.files("./civet_run", full.names = T)

#Dimension names for array
sample_names <- rownames(best.genoprobs.192)
All_transcriptID_list <- read.delim(file = paste(civet_dir[1],"/gbrs.quantified.diploid.isoforms.tpm", sep =""), header = TRUE, sep ="\t")
All_transcriptID_list <- All_transcriptID_list$locus

#Create a data frame with all transcripts 
All_transcript_tpm <- array(0, c(length(sample_names), length(All_transcriptID_list)), dimnames = list(sample_names, All_transcriptID_list))
for (i in 1:length(sample_names)){
	temp <- read.delim(file = paste(civet_dir[i],"/gbrs.quantified.diploid.isoforms.tpm", sep =""), header = TRUE, sep ="\t")
	rownames(temp) <- temp$locus
	All_transcript_tpm[i, ] <- temp$total 
}
save(All_transcript_tpm, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")

#Get RankZ score
All_transcript_tpm <- t(All_transcript_tpm)
All_transcript_tpm_rankZ <- apply(All_transcript_tpm_rankZ, 2, rankZ)
All_transcript_tpm_rankZ <- t(All_transcript_tpm_rankZ)
save(All_transcript_tpm_rankZ, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")

#Extract Fmn1 transcripts of all isoforms
Fmn1_transcriptID_list <- c("ENSMUST00000102547", "ENSMUST00000081349", "ENSMUST00000110954", "ENSMUST00000152255", "ENSMUST00000153151", "ENSMUST00000145891",
					"ENSMUST00000150510", "ENSMUST00000161731", "ENSMUST00000099576")
#Create new dataframe with just Fmn1 RankZ counts

Fmn1_transcript_tpm_rankZ <- All_transcript_tpm_rankZ[,Fmn1_transcriptID_list]
save(Fmn1_transcript_tpm_rankZ, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Fmn1_transcript_tpm_rankZ.Rdata")
Fmn1_transcript_tpm <- All_transcript_tpm[,Fmn1_transcriptID_list]
save(Fmn1_transcript_tpm, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Fmn1_transcript_tpm.Rdata")

###############################################
#Plot data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Fmn1_transcript_tpm_rankZ.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Fmn1_transcript_tpm.Rdata")
Fmn1_names <- c("Fmn1_001", "Fmn1_002", "Fmn1_003", "Fmn1_004", "Fmn1_005", "Fmn1_006", "Fmn1_007", "Fmn1_008", "Fmn1_201")

#plot dot plot with discrete x and continuous y
#change dataframe format 
gg_data <- Fmn1_transcript_tpm_rankZ
colnames(gg_data) <- Fmn1_names
gg_data <- melt(gg_data, id.vars = Fmn1_names,
	value.name = "tpm_rankZ")
names(gg_data)[2] <- "Fmn1_transcripts"

#plot RankZ tpm
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcirpt_tpm_rankZ.png", width = 1500, height = 1000, res = 100)
ggplot(gg_data, aes(x =Fmn1_transcripts, y = tpm_rankZ, fill =Fmn1_transcripts)) + 
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.01) +
	scale_x_discrete("Fmn1 transcirpts") + 
	scale_y_continuous("transcript tpm rankZ") +
	labs( title = "Comparison of Fmn1 transcirpts tpm rankZ") +
	theme( plot.title = element_text(hjust = 0.5))
dev.off()

#plot RankZ tpm by allele
gg_data <- Fmn1_transcript_tpm_rankZ
gg_data <- as.data.frame(gg_data)
colnames(gg_data) <- Fmn1_names
gg_data$allele <- Gene_allele$ENSMUSG00000044042
gg_data <- melt(gg_data)
colnames(gg_data) <- c("Allele", "Fmn1_transcripts", "Value")

ggplot <- ggplot(gg_data, aes( x = Fmn1_transcripts, y = Value, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.03, position = position_dodge(1.5)) +
	scale_x_discrete("Fmn1 transcirpts") + 
	scale_y_continuous("transcript rankZ tpm") +
	labs( title = "Comparison of Fmn1 transcirpts rankZ tpm by allele") +
	theme( plot.title = element_text(hjust = 0.5))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcript_rankZ_tpm_by_allele.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

#plot non-transformed tpm
gg_data <- Fmn1_transcript_tpm
gg_data <- as.data.frame(gg_data)
colnames(gg_data) <- Fmn1_names
gg_data <- melt(gg_data, id.vars = Fmn1_names,
	value.name = "tpm")
names(gg_data)[2] <- "Fmn1_transcripts"

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcirpt_tpm.png", width = 1500, height = 1000, res = 100)
ggplot(gg_data, aes(x =Fmn1_transcripts, y = tpm, fill =Fmn1_transcripts)) + 
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.1) +
	scale_x_discrete("Fmn1 transcirpts") + 
	scale_y_continuous("transcript tpm") +
	labs( title = "Comparison of Fmn1 transcirpts tpm") +
	theme( plot.title = element_text(hjust = 0.5))
dev.off()

#plot non-transformed tpm by allele 
gg_data <- Fmn1_transcript_tpm
gg_data <- as.data.frame(gg_data)
colnames(gg_data) <- Fmn1_names
gg_data$allele <- Gene_allele$ENSMUSG00000044042
gg_data <- melt(gg_data)
colnames(gg_data) <- c("Allele", "Fmn1_transcripts", "Value")

ggplot <- ggplot(gg_data, aes( x = Fmn1_transcripts, y = Value, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2, position = position_dodge(1)) +
	scale_x_discrete("Fmn1 transcirpts") + 
	scale_y_continuous("transcript tpm") +
	labs( title = "Comparison of Fmn1 transcirpts tpm by allele") +
	theme( plot.title = element_text(hjust = 0.5))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

#plot non-transformed tpm by allele but remove non-coding transcripts
gg_data <- Fmn1_transcript_tpm
gg_data <- as.data.frame(gg_data)
colnames(gg_data) <- Fmn1_names
gg_data$allele <- Gene_allele$ENSMUSG00000044042
gg_data$Fmn1_003 <- NULL
gg_data$Fmn1_004 <- NULL
gg_data$Fmn1_005 <- NULL
gg_data$Fmn1_006 <- NULL
gg_data$Fmn1_007 <- NULL
gg_data <- melt(gg_data)
colnames(gg_data) <- c("Allele", "Fmn1_transcripts", "Value")

ggplot <- ggplot(gg_data, aes( x = Fmn1_transcripts, y = Value, fill = Fmn1_transcripts)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2) +
	scale_x_discrete("Fmn1 transcirpts") + 
	scale_y_continuous("transcript tpm") +
	labs( title = "Comparison of Fmn1 transcirpts tpm by allele") +
	theme( plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_coding_transcript_tpm.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

ggplot <- ggplot(gg_data, aes( x = Fmn1_transcripts, y = Value, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2, position = position_dodge(1)) +
	scale_x_discrete("Fmn1 transcirpts") + 
	scale_y_continuous("transcript tpm") +
	labs( title = "Comparison of Fmn1 transcirpts tpm by allele") +
	theme( plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_coding_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

#eQTL of 3 transcripts with most expression
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/sex.covar.Rdata")

#Fmn1-002
#qtl.Fmn1_002 <- scanone( pheno = All_transcript_tpm_rankZ, pheno.col = "ENSMUST00000081349", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#save(qtl.Fmn1_002, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1_002.Rdata")
#Fmn1-005
qtl.Fmn1_005 <- scanone( pheno = All_transcript_tpm_rankZ, pheno.col = "ENSMUST00000153151", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.Fmn1_005, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1_005.Rdata")
#Fmn1-201
qtl.Fmn1_201 <- scanone( pheno = All_transcript_tpm_rankZ, pheno.col = "ENSMUST00000099576", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
save(qtl.Fmn1_201, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1_201.Rdata")

#plot eQTL
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1_002.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1_005.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1_201.Rdata")

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_002.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Fmn1_002, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Fmn1_002 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_005.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Fmn1_005, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Fmn1_005 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_201.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Fmn1_201, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Fmn1_201 QTL perm 1000")
dev.off()

#	Allele effect plots
qtl.Fmn1_002$coef$A[abs(qtl.Fmn1_002$coef$A) > 0.05 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_002.coef.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(qtl.Fmn1_002, chr = 2, main = "Allele effect of Fmn1_002 QTL at Chr 2")
dev.off()
qtl.Fmn1_005$coef$A[abs(qtl.Fmn1_005$coef$A) > 0.05 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_005.coef.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(qtl.Fmn1_005, chr = 2, main = "Allele effect of Fmn1_005 QTL at Chr 2")
dev.off()
qtl.Fmn1_201$coef$A[abs(qtl.Fmn1_201$coef$A) > 0.05 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1_201.coef.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(qtl.Fmn1_201, chr = 2, main = "Allele effect of Fmn1_201 QTL at Chr 2")
dev.off()







