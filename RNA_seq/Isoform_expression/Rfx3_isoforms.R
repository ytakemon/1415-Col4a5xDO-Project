#Col4a5xDO Rfx isoforms
#Yuka Takemon
#Created: 12/09/16
#Last modified: 12/09/16

#Gene Ensembl ID: ENSMUSG00000040929

#Transcripts according to Ensembl:
#Rfx3_001: ENSMUST00000174850
#Rfx3_002: ENSMUST00000173161
#Rfx3_003: ENSMUST00000046898 
#Rfx3_004: ENSMUST00000172907
#Rfx3_005: ENSMUST00000173863
#Rfx3_006: ENSMUST00000172498
#Rfx3_007: ENSMUST00000174420
#Rfx3_201: ENSMUST00000165566


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



