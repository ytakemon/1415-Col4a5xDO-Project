################################################################################################################################
## Col4a5xDO Tgm2 isoform sequence
## Yuka Takemon
## Created: 08/03/17

## OBJECTIVE


## sessionInfo()

## List of data saved from this script (time savers for reanalysis)
## Rdata:

## Tables:
## Plots:

## Ensembl ID for genes:
# Tgm2: ENSMUSG00000037820

## Tgm2 transcript IDs:
# Tgm2-201 ENSMUST00000103122
# Tgm2-202 ENSMUST00000140923
# Tgm2-204 ENSMUST00000174718
# Tgm2-203 ENSMUST00000152690
################################################################################################################################
library(DOQTL)
library(ggplot2)
library(reshape2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load essential data and tailor to Tgm2 -----------------------------------------
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
Gene_allele <- as.data.frame(Gene_allele)
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")

#Extract Tgm2 transcripts of all isoforms
Tgm2_transcriptID_list <- c("ENSMUST00000103122", "ENSMUST00000140923", "ENSMUST00000174718", "ENSMUST00000152690")

#Create new dataframe with just Tgm2 RankZ counts
Tgm2_transcript_tpm_rankZ <- All_transcript_tpm_rankZ[,Tgm2_transcriptID_list]
save(Tgm2_transcript_tpm_rankZ, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Tgm2_transcript_tpm_rankZ.Rdata")
Tgm2_transcript_tpm <- All_transcript_tpm[,Tgm2_transcriptID_list]
save(Tgm2_transcript_tpm, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Tgm2_transcript_tpm.Rdata")

#Plot data ---------------------------------------------------------------------
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Tgm2_transcript_tpm_rankZ.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Tgm2_transcript_tpm.Rdata")
Tgm2_names <- c("ENSMUST00000103122", "ENSMUST00000140923", "ENSMUST00000174718", "ENSMUST00000152690")

#plot dot plot with discrete x and continuous y
#change dataframe format
gg_data <- Tgm2_transcript_tpm_rankZ
colnames(gg_data) <- Tgm2_names
gg_data <- melt(gg_data, id.vars = Tgm2_names,
	value.name = "tpm_rankZ")
names(gg_data)[2] <- "Tgm2_transcripts"

#plot RankZ tpm
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tgm2_transcirpt_tpm_rankZ.pdf", width = 7, height = 5)
ggplot(gg_data, aes(x =Tgm2_transcripts, y = tpm_rankZ, fill =Tgm2_transcripts)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.02) +
	scale_x_discrete("Tgm2 transcirpts") +
	scale_y_continuous("Rankz Transformed TPM") +
	labs( title = "Comparison of Tgm2 isoform TPM") +
  guides(fill = FALSE)+
	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
dev.off()

#plot RankZ tpm by allele
gg_data <- Tgm2_transcript_tpm_rankZ
gg_data <- as.data.frame(gg_data)
colnames(gg_data) <- Tgm2_names
gg_data$allele <- Gene_allele$ENSMUSG00000037820
gg_data <- melt(gg_data)
colnames(gg_data) <- c("Allele", "Tgm2_transcripts", "Value")

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tgm2_transcript_rankZ_tpm_by_allele.pdf", width = 10.0, height = 7.5)
ggplot(gg_data, aes( x = Tgm2_transcripts, y = Value, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.03, position = position_dodge(1.5)) +
	scale_x_discrete("Tgm2 transcirpts") +
	scale_y_continuous("RankZ transformed TPM") +
	labs( title = "Comparison of Tgm2 transcirpts rankZ TPM by allele") +
	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
dev.off()

#plot non-transformed tpm
gg_data <- Tgm2_transcript_tpm
gg_data <- melt(gg_data, id.vars = Tgm2_names,
	value.name = "tpm")
names(gg_data)[2] <- "Tgm2_transcripts"

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tgm2_transcirpt_tpm.pdf", width = 7, height = 5)
ggplot(gg_data, aes(x =Tgm2_transcripts, y = tpm, fill =Tgm2_transcripts)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 3) +
	scale_x_discrete("Tgm2 transcirpts") +
	scale_y_continuous("Transcript TPM") +
  guides(fill = FALSE)+
	labs( title = "Comparison of Tgm2 transcript TPM") +
	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
dev.off()

#Plot non-transformd tpm by allele
gg_data <- Tgm2_transcript_tpm
gg_data <- as.data.frame(gg_data)
colnames(gg_data) <- Tgm2_names
gg_data$allele <- Gene_allele$ENSMUSG00000037820
gg_data <- melt(gg_data)
colnames(gg_data) <- c("Allele", "Tgm2_transcripts", "Value")

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tgm2_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)
ggplot(gg_data, aes( x = Tgm2_transcripts, y = Value, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 3, position = position_dodge(1.5)) +
	scale_x_discrete("Tgm2 transcirpts") +
	scale_y_continuous("Transcript TPM") +
	labs( title = "Comparison of Tgm2 transcirpts TPM by allele") +
	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
dev.off()






# Make allele effect plots of Tgm2 eQTL ----------------------------------------
#eQTL of 3 transcripts with most expression -
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm_rankZ.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/sex.covar.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/ENSMUSG00000037820.Tgm2.eQTL.rankZ.tpm.Rdata")
qtl <- ENSMUSG00000037820.Tgm2.eQTL

# Plot whole chromosome 2
qtl$coef$A[abs(qtl$coef$A) > 0.1 ] = 0
coefplot(qtl, chr = 2, main = "Allele effect of Tgm2_002 QTL at Chr 2")

#	subset qtl object to region of interest
# LOD
lodA <- qtl$lod$A
lodA <- lodA[lodA$chr == 2,]
lodA <- lodA[lodA$pos > 140, ]
lodA <- lodA[lodA$pos < 170, ]
# Coef
coefA <- qtl$coef$A #matrix
coefA <- as.data.frame(coefA)
coefA <- coefA[rownames(lodA),]
# replace
qtl$lod$A <- lodA
qtl$coef$A <- coefA
qtl$coef$A[abs(qtl$coef$A) > 0.1 ] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tgm2_allele_effect_chr2.pdf", width = 7, height = 5)
coefplot(qtl, chr = 2, main = "Allele effect of Tgm2 QTL at Chr 2")
dev.off()



## Plotting narrower allele effect Alb 10wk
#	subset qtl object to region of interest
# LOD
lodA <- qtl$lod$A
lodA <- lodA[lodA$chr == 2,]
lodA <- lodA[lodA$pos > 101.7, ]
lodA <- lodA[lodA$pos < 113.7, ]
# Coef
coefA <- qtl$coef$A #matrix
coefA <- as.data.frame(coefA)
coefA <- coefA[rownames(lodA),]
# replace
qtl$lod$A <- lodA
qtl$coef$A <- coefA
#	plot allele effect
qtl$coef$A[abs(qtl$coef$A) > 2 ] = 0
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr2.bayesian.int.pdf", width = 10.0, height = 7.5)
coefplot(qtl, chr = 2, main = "Chr 2 Alb10WK_log Allele Effect Plot @ bayesian interval")
dev.off()
