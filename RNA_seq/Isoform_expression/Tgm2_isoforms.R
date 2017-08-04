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

#ANOVA of Tgm2 TPM by allele
gg_data <- Tgm2_transcript_tpm
gg_data <- cbind(gg_data, as.character(Gene_allele$ENSMUSG00000037820))
colnames(gg_data)[5] <- "Allele"
gg_melt <- gg_data[,1:4]
gg_melt <- melt(gg_melt, id.vars = Tgm2_names,
					value.name = "tpm")
names(gg_melt)[2] <- "Tgm2_transcripts"
Tgm2_names4 <- c(gg_data[,5], gg_data[,5], gg_data[,5], gg_data[,5])
gg_melt <- cbind(gg_melt, Tgm2_names4)
gg_melt <- as.data.frame(gg_melt, stringsAsFactors = FALSE)
gg_melt <- gg_melt[gg_melt$Tgm2_transcripts %in% "ENSMUST00000103122",]
gg_melt$tpm <- as.numeric(as.character(gg_melt$tpm))

fit <- lm(tpm ~ Tgm2_names4, data = gg_melt)
summary(fit)

#Call:
#lm(formula = tpm ~ Tgm2_names4, data = gg_melt)
#
#Residuals:
#   Min     1Q Median     3Q    Max
#-66.41 -26.32 -10.40  22.72 129.39
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)   134.7785     7.9592  16.934   <2e-16 ***
#Tgm2_names4BB  -6.6057    10.6359  -0.621   0.5353
#Tgm2_names4BC   0.8770    10.3393   0.085   0.9325
#Tgm2_names4BD  -0.7576    10.2146  -0.074   0.9410
#Tgm2_names4BE  19.5405    11.6918   1.671   0.0964 .
#Tgm2_names4BF   7.2041    12.5004   0.576   0.5651
#Tgm2_names4BG -34.8265    15.4129  -2.260   0.0250 *
#Tgm2_names4BH -20.1577    10.2146  -1.973   0.0499 *
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Residual standard error: 37.33 on 184 degrees of freedom
#Multiple R-squared:  0.1063,	Adjusted R-squared:  0.07229
#F-statistic: 3.126 on 7 and 184 DF,  p-value: 0.003839
colnames(gg_melt) <- c("SampleID", "Transcripts", "TPM", "Allele")
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tgm2_transcript_boxplot.pdf", width = 7, height = 7.5)
ggplot( gg_melt, aes( x = Allele, y = TPM, fill = Allele)) +
      geom_boxplot() +
      scale_x_discrete() +
      theme( text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
      xlab( "Allele") +
      ylab( "TPM") +
      ggtitle( "ENSMUST00000103122 by DO diet groups")
dev.off()

#Reference AB
gg_melt_AB <- within(gg_melt, Allele <- relevel(Allele, ref = "AB"))
summary(lm(TPM ~ Allele, data = gg_melt_AB))
# sig anova, pair: BG, BH
#Reference BB
gg_melt_BB <- within(gg_melt, Allele <- relevel(Allele, ref = "BB"))
summary(lm(TPM ~ Allele, data = gg_melt_BB))
# sig anova, pair: BE
#Reference BC
gg_melt_BC <- within(gg_melt, Allele <- relevel(Allele, ref = "BC"))
summary(lm(TPM ~ Allele, data = gg_melt_BC))
# sig anova, pair: BG, BH
#Reference BD
gg_melt_BD <- within(gg_melt, Allele <- relevel(Allele, ref = "BD"))
summary(lm(TPM ~ Allele, data = gg_melt_BD))
# sig anova, pair: BG,BH
#Reference BE
gg_melt_BE <- within(gg_melt, Allele <- relevel(Allele, ref = "BE"))
summary(lm(TPM ~ Allele, data = gg_melt_BE))
# sig anova, pair: BB, BG, BH
#Reference BF
gg_melt_BF <- within(gg_melt, Allele <- relevel(Allele, ref = "BF"))
summary(lm(TPM ~ Allele, data = gg_melt_BF))
# sig anova, pair: BG, BH
#Reference BG
gg_melt_BG <- within(gg_melt, Allele <- relevel(Allele, ref = "BG"))
summary(lm(TPM ~ Allele, data = gg_melt_BG))
# sig anova, pair: AB, BE, BF
#Reference BH
gg_melt_BH <- within(gg_melt, Allele <- relevel(Allele, ref = "BH"))
summary(lm(TPM ~ Allele, data = gg_melt_BH))
# sig anova, pair: AB, BC, BD, BE, BF



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
