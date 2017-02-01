#Create plot of Rictor expression with GFR, Alb by allele
# phenotype on y-axis, Rictor expression on x-axis
# add colours by Gene allele
# Calculate average of GFR exp by Rictor allele and do ANOVA

#Rictor ENSMUSG00000050310 38exons chr15

# load library and data
library(DOQTL)
library(ggplot2)
library(reshape2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clearn up pheno data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
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
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Creat6WK_log","Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log", "ACR6WK_log", "ACR10WK_log", "ACR15WK_log")]

#clearn up pheno data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192

#add Rcitor into pheno ENSMUSG00000050310 using RankZ-transformed tpm values
RNA_seqZ <- as.data.frame(RNA_seqZ)
Gene_allele <- as.data.frame(Gene_allele)
pheno$Rictor_tpmZ <- RNA_seqZ$ENSMUSG00000050310
pheno$Rictor_allele <- Gene_allele$ENSMUSG00000050310
pheno_M <- pheno[pheno$Sex == "M",]
options(na.action = 'na.pass')

#Correlation between Rictor x GFR
G_pheno <- pheno[complete.cases(pheno$C2_log),]
fit <- lm(C2_log ~ Rictor_tpmZ , G_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(G_pheno, aes(y = C2_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 7.5, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed C2 GFR vs RankZ-transformed Rictor TPM") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvC2GFR.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off() 
G_pheno <- pheno_M[complete.cases(pheno_M$C2_log),]
fit <- lm(C2_log ~ Rictor_tpmZ , G_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(G_pheno, aes(y = C2_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 7.5, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed C2 GFR vs RankZ-transformed Rictor TPM (Males)") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvC2GFR_Males.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

#Correlation between Rictor x ACR6
A6_pheno <- pheno[complete.cases(pheno$ACR6WK_log),]
fit <- lm(ACR6WK_log ~ Rictor_tpmZ, A6_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(A6_pheno, aes(y = ACR6WK_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 8.8, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed ACR6WK", breaks = seq(0, 8.8, by = 0.5)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed ACR6WK vs RankZ-transformed Rictor TPM") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvACR6WK.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()
A6_pheno <- pheno_M[complete.cases(pheno_M$ACR6WK_log),]
fit <- lm(ACR6WK_log ~ Rictor_tpmZ, A6_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(A6_pheno, aes(y = ACR6WK_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 8.8, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed ACR6WK", breaks = seq(0, 8.8, by = 0.5)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed ACR6WK vs RankZ-transformed Rictor TPM (Males)") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvACR6WK_males.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

#Correlation between Rictor x ACR10
A10_pheno <- pheno[complete.cases(pheno$ACR10WK_log),]
fit <- lm(ACR10WK_log ~ Rictor_tpmZ, A10_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(A10_pheno, aes(y = ACR10WK_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 8.8, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed ACR10WK", breaks = seq(0, 8.8, by = 0.5)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed ACR10WK vs RankZ-transformed Rictor TPM") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvACR10WK.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()
A10_pheno <- pheno_M[complete.cases(pheno_M$ACR10WK_log),]
fit <- lm(ACR10WK_log ~ Rictor_tpmZ, A10_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(A10_pheno, aes(y = ACR10WK_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 8.8, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed ACR10WK", breaks = seq(0, 8.8, by = 0.5)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed ACR10WK vs RankZ-transformed Rictor TPM (Males)") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvACR10WK_males.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

#Correlation between Rictor x ACR15
A15_pheno <- pheno[complete.cases(pheno$ACR15WK_log),]
fit <- lm(ACR15WK_log ~ Rictor_tpmZ, A15_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(A15_pheno, aes(y = ACR15WK_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 9.8, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed ACR15WK", breaks = seq(0, 9.8, by = 0.5)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed ACR15WK vs RankZ-transformed Rictor TPM") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvACR15WK.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()
A15_pheno <- pheno_M[complete.cases(pheno_M$ACR15WK_log),]
fit <- lm(ACR15WK_log ~ Rictor_tpmZ, A15_pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(A15_pheno, aes(y = ACR15WK_log, x = Rictor_tpmZ))+
	geom_point(size = 2, aes(colour = Rictor_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 9.8, x = 0.925, label = eq, fontface = "bold", size = 5) + 
	scale_y_continuous( "Log-transformed ACR15WK", breaks = seq(0, 9.8, by = 0.5)) +
	scale_x_continuous( "RankZ-transformed Rictor TPM") +
	labs( title = "Log-transformed ACR15WK vs RankZ-transformed Rictor TPM (Males)") +
	theme(plot.title = element_text(hjust = 0.5))  
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_RictorvACR15WK_males.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

# Calculate average of GFR exp by Rictor allele and do ANOVA
#plots GFR by founder
G_pheno <- pheno_M[complete.cases(pheno_M$C2_log),]
ggplot <- ggplot(G_pheno, aes(x = Rictor_allele, y = C2_log)) +
	geom_boxplot( fill = "grey80", colour = "black") +
	scale_x_discrete("Rictor DO founder allele") +
	scale_y_continuous("Log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) + 
	labs( title = "Log-transformed C2 GFR by DO founder (Males)")+
	theme(plot.title = element_text(hjust = 0.5))
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GFRbyFounder_males.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()
#calcualte anova
library(knitr)
GFR_aov <- aov(C2_log ~ Rictor_allele, data = G_pheno)
knitr::kable((anova(GFR_aov)))
#|              | Df|   Sum Sq|   Mean Sq|   F value|    Pr(>F)|
#|:-------------|--:|--------:|---------:|---------:|---------:|
#|Rictor_allele |  6|  1.78789| 0.2979817| 0.5503983| 0.7674729|
#|Residuals     | 54| 29.23521| 0.5413928|        NA|        NA|

#no significant differences, but ron said it mimics trends in human cyprus data.










