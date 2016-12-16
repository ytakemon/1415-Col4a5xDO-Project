#Plot corelation between Rfx3 expression and phenotype
#Yuka Takemon
#12/02/2016

#Ensembl ID for genes:
#Rfx3: ENSMUSG00000040929 ENSMUSG00000040929
#Fmn1: ENSMUSG00000044042

library(DOQTL)
library(ggplot2)
library(grid)

setwd("/hpcdata/ytakemon/Col4a5xDO")
load("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(RNA_seq),] #subset pheno to match 192 samples
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$C2_log <- log(pheno$C2) 
pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR10WK_log <- log(pheno$ACR10WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)
pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Alb10WK_log <- log(pheno$Alb10WK)
pheno$Alb15WK_log <- log(pheno$Alb15WK)
pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno$Creat10WK_log <- log(pheno$Creat10WK)
pheno$Creat15WK_log <- log(pheno$Creat15WK)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
#keep only interested columns
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Creat6WK_log",
	"Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log",
	"ACR6WK_log","ACR10WK_log","ACR15WK_log")]

#subset out Rfx3 and Fmn1
RNA_subset <- RNA_seq[, c("ENSMUSG00000040929", "ENSMUSG00000044042")]
colnames(RNA_subset) <- c("Rfx3", "Fmn1")
RNA_subset <- as.data.frame(RNA_subset)
#consolidate RNA with pheno
RNA_pheno <- pheno
RNA_pheno$Rfx3 <- RNA_subset$Rfx3
RNA_pheno$Fmn1 <- RNA_subset$Fmn1

#calculate and plot correlation
#Rfx3
data <- RNA_pheno[ complete.cases(RNA_pheno$C2_log),]

fit <- lm(Rfx3 ~ C2_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
eq2 <- paste(pval)
gfr <- ggplot(data, aes(y = Rfx3, x = C2_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(size = 2)) +
	annotate( "text" , y = 2.6, x = 4.6, label = eq, fontface = "bold", size = 3) + 
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Rfx3 tpm", breaks = seq(0, 7.5, by = 0.2)) +
	labs( title = "Rfx3 TPM vs log C2 GFR") +
	theme(plot.title = element_text(hjust = 0.5))  

data <- RNA_pheno[ complete.cases(RNA_pheno$ACR6WK_log),]

fit <- lm(Fmn1 ~ ACR6WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
eq2 <- paste(pval)
ACR <- ggplot(data, aes(y = Rfx3, x = ACR6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(size = 2)) +
	annotate( "text" , y = 2.6, x = 4.6, label = eq, fontface = "bold", size = 3) + 
	scale_x_continuous( "log-transformed ACR6wk", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Rfx3 tpm", breaks = seq(0, 7.5, by = 0.2)) +
	labs( title = "Rfx3 TPM vs log ACR6wk") +
	theme(plot.title = element_text(hjust = 0.5))  

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/test.png", width = 2000, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(1, 2)))
print(gfr, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()


#need to add allele to colour in dots for Rfx3 and Fmn1



