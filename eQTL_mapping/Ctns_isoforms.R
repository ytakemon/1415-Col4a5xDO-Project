#Col4a5xDO Ctns isoforms
#Yuka Takemon
#Created: 12/09/16
#Last modified: 12/09/16

#Gene Ensembl ID:
#Ctns: ENSMUSG00000005949

#Transcripts according to Ensembl:
#Ctns_001: ENSMUST00000006103.8 (11 exons, protein coding)
#Ctns_002: ENSMUST00000108476.7 (12 exons, protein coding)
#Ctns_003: ENSMUST00000150468.1 (6 exons, processed transcript) 
#Ctns_004: ENSMUST00000144658.7 (4 exons, retained intron)
#Ctns_005: ENSMUST00000130101.1 (4 exons, retained intron)

#Transcript ID as found in isoform.tpm files
#Ctns_001: ENSMUST00000006103
#Ctns_002: ENSMUST00000108476
#Ctns_003: ENSMUST00000150468
#Ctns_004: ENSMUST00000144658
#Ctns_005: ENSMUST00000130101

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
Ctns_ID <- c("ENSMUST00000006103", "ENSMUST00000108476", "ENSMUST00000150468", "ENSMUST00000144658", "ENSMUST00000130101")
Ctns_names <- c("Ctns_001", "Ctns_002", "Ctns_003", "Ctns_004", "Ctns_005")
Gene_allele <- as.data.frame(Gene_allele)

All_transcript_tpm <- as.data.frame(All_transcript_tpm)
All_transcript_tpm_rankZ <- as.data.frame(All_transcript_tpm_rankZ)

Ctns_ID_tpm <- All_transcript_tpm[, Ctns_ID]
Ctns_ID_ranked_tpm <- All_transcript_tpm_rankZ[, Ctns_ID]

gg_data <- Ctns_ID_tpm
colnames(gg_data) <- Ctns_names
gg_data$allele <- Gene_allele$ENSMUSG00000005949
gg_data <- melt(gg_data)
names(gg_data) <- c("Allele", "Ctns_transcript_ID", "Values")

ggplot <- ggplot(gg_data, aes( x = Ctns_transcript_ID, y = Values, fill = Ctns_transcript_ID)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2) +
	scale_x_discrete( "Ctns transcripts") +
	scale_y_continuous( "transcript tpm") +
	labs( title = "Comparison of Ctns transcirpts tpm") +
	theme( plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Ctns_transcript_tpm.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

ggplot <- ggplot(gg_data, aes( x = Ctns_transcript_ID, y = Values, fill = Allele)) +
	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2, position = position_dodge(1)) +
	scale_x_discrete( "Ctns transcripts") +
	scale_y_continuous( "transcript tpm") +
	labs( title = "Comparison of Ctns transcirpts tpm") +
	theme( plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Ctns_transcript_tpm_by_allele.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()



