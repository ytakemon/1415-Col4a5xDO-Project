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


######
fit <- lm(Col4a2_tpm ~ Col4a1_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot12 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a2_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , x = 1.9, y = 2.3, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a1 vs Col4a2 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 