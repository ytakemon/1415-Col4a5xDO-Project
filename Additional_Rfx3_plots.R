##########################################################################################################################################################
## Additional plots for Rfx3																															##
## Contributor: Yuka Takemon 																															##
## Created: 01/24/17																																	##
##																																						##
## OBJECTIVE:																																			##
##	Allele effect plot of Rfx3 expression at Rfx3 locus 																								##
##	Allele effect plot of GFR at Rfx3 locus 																											##
##	Correlation between Rfx3 exp vs .....Dync2li1, Foxj1, Dnahc5, Dnahc11, RankZ_Rfx3vDnahc9 															##
##	Coefficeint plot of Rfx3 narrower focus on LOD score peak																							##
## 																																						##
## RESULTS:																																				##
##																																						##
##########################################################################################################################################################

###############################################################################################################################
## Additional plots for Rfx3
#	Name of file will probably change 
#	Allele effect plot of Rfx3 expression at Rfx3 locus
#	Allele effect plot of GFR at Rfx3 locus
#	Correlation between Rfx3 exp vs .....Dync2li1, Foxj1, Dnahc5, Dnahc11, Dnahc9
#	Coefficeint plot of Rfx3 narrower focus on LOD score peak
## Yuka Takemon
## 01/24/17

# load library and data
library(DOQTL)
library(ggplot2)
library(knitr)
library(reshape2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Rfx3.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clearn up pheno data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2) 
#add Rfx3 into pheno ENSMUSG00000040929 (using raw tpm reads for allele effect)
RNA_seq <- as.data.frame(RNA_seq)
pheno$Rfx3_tpm <- RNA_seq$ENSMUSG00000040929
options(na.action = 'na.pass')

###############################################################################################################

## Allele effect plot of Rfx3 expression at Rfx3 locus#######################################################
#	find snp that falls on Rfx3
data <- GM_snps[GM_snps$chr == 19,]
data <- data[data$pos > 27.70969, ]
data <- data[data$pos < 28.10969, ]
#chosen snp: UNC30148273, rs6279232

fit <- lm(pheno$Rfx3_tpm ~ pheno$Sex + best.genoprobs.192[,, "UNC30148273"], na.action = na.exclude)
#Female averages by strain at Snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1]  #add female 
#ecfit <- exp(cfit) #exp fixes values back to natural numbers
cfit_F <- as.data.frame(cfit)
#Male averages by strain at snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1] + cfit[2] #add female and male = male (female = female, male = male + female)
#ecfit <- exp(cfit) #exp fixes balues back to natural numbers
cfit_M <- as.data.frame(cfit)

data_coef <- cfit_F
data_coef$M <- cfit_M$cfit
colnames(data_coef) <- c("F","M")
names(data_coef) <- c("F","M")
data_coef <- data_coef[c(3:10),]
rownames(data_coef) <- c("A","B","C","D","E","F","G","H")
data_coef$founder <- rownames(data_coef)

ggdata <- melt(data_coef)
colnames(ggdata) <- c("Founder", "Sex", "Value")
ggplot <- ggplot(ggdata, aes( x = Founder, y = Value, fill = Sex)) +
	geom_bar( stat = "identity", position = "dodge") +
	scale_fill_discrete( labels = c("Females","Males")) +
	labs(title = "Col4a5xDO Rfx3 (ENSMUSG00000040929) TPM @  Rfx3 locus UNC30148273 (rs6279232)", x = "DO Founders", y = "Rfx3 TPM") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3_TPM_at_Rfx3_UNC30148273_rs6279232_by_Founders.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

## Allele effect plot of GFR at Rfx3 locus#####################################################################
fit <- lm(pheno$C2_log ~ pheno$Sex + best.genoprobs.192[,, "UNC30148273"], na.action = na.exclude)
#Female averages by strain at Snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1]  #add female 
#ecfit <- exp(cfit) #exp fixes values back to natural numbers
cfit_F <- as.data.frame(cfit)
#Male averages by strain at snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1] + cfit[2] #add female and male = male (female = female, male = male + female)
#ecfit <- exp(cfit) #exp fixes balues back to natural numbers
cfit_M <- as.data.frame(cfit)

data_coef <- cfit_F
data_coef$M <- cfit_M$cfit
colnames(data_coef) <- c("F","M")
names(data_coef) <- c("F","M")
data_coef <- data_coef[c(3:10),]
rownames(data_coef) <- c("A","B","C","D","E","F","G","H")
data_coef$founder <- rownames(data_coef)

ggdata <- melt(data_coef)
colnames(ggdata) <- c("Founder", "Sex", "Value")
ggplot <- ggplot(ggdata, aes( x = Founder, y = Value, fill = Sex)) +
	geom_bar( stat = "identity", position = "dodge") +
	scale_fill_discrete( labels = c("Females","Males")) +
	labs(title = "Col4a5xDO log-transformed C2 GFR @ Rfx3 locus UNC30148273 (rs6279232)", x = "DO Founders", y = "Log C2 GFR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/C2_GFR_at_Rfx3_UNC30148273_rs6279232_by_Founders.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()


## Correlation between Rfx3 exp vs .....Dync2li1, Foxj1, Dnahc5, Dnahc11, Dnahc9 #########################################################
#Rfx3 ENSMUSG00000040929
#Dync2li1 ENSMUSG00000024253
#Foxj1 ENSMUSG00000034227
#Dnahc5 ENSMUSG00000022262
#Dnahc11 ENSMUSG00000018581
#Dnahc9 ENSMUSG00000056752

#Use rankz RNA-seq for this because we are comparing between rna-seq tpms
RNA_seqZ <- as.data.frame(RNA_seqZ)

pheno$Rfx3_ranked <- RNA_seqZ$ENSMUSG00000040929
pheno$Dync2li1 <- RNA_seqZ$ENSMUSG00000024253
pheno$Foxj1 <- RNA_seqZ$ENSMUSG00000034227
pheno$Dnahc5 <- RNA_seqZ$ENSMUSG00000022262
pheno$Dnahc11 <- RNA_seqZ$ENSMUSG00000018581
pheno$Dnahc9 <- RNA_seqZ$ENSMUSG00000056752

# Rfx3 v. Dync2li1
fit <- lm(Rfx3 ~ Dync2li1, pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno, aes(y = Rfx3, x = Dync2li1)) +
	geom_point(size = 2) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 0.7, x = 1.65, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-tranformed Dync2li1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 0.7, by = 0.01)) +
	labs( title = "Rfx3 TPM vs Dync2li1 TPM") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_Rfx3vDync2li1.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

# Rfx3 v. Foxj1
fit <- lm(Rfx3 ~ Foxj1, pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno, aes(y = Rfx3, x = Foxj1)) +
	geom_point(size = 2) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 0.7, x = 0.6, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-tranformed Foxj1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 0.7, by = 0.01)) +
	labs( title = "Rfx3 TPM vs Foxj1 TPM") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_Rfx3vFoxj1.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

# Rfx3 v. Dnahc5
fit <- lm(Rfx3 ~ Dnahc5, pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno, aes(y = Rfx3, x = Dnahc5)) +
	geom_point(size = 2) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 0.7, x = 0.5, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-tranformed Dnahc5 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 0.7, by = 0.01)) +
	labs( title = "Rfx3 TPM vs Dnahc5 TPM") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_Rfx3vDnahc5.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

# Rfx3 v. Dnahc11
fit <- lm(Rfx3 ~ Dnahc11, pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno, aes(y = Rfx3, x = Dnahc11)) +
	geom_point(size = 2) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 0.7, x = 0.5, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-tranformed Dnahc11 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 0.7, by = 0.01)) +
	labs( title = "Rfx3 TPM vs Dnahc11 TPM") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_Rfx3vDnahc11.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

# Rfx3 v. Dnahc9
# There is a negative value because TPM pre transformation was 0
# correct tom back to 0
pheno$Dnahc9[pheno$Dnahc9 < 0] = 0

fit <- lm(Rfx3 ~ Dnahc9, pheno)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno, aes(y = Rfx3, x = Dnahc9)) +
	geom_point(size = 2) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 0.7, x = 0.3, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-tranformed Dnahc9 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 0.7, by = 0.01)) +
	labs( title = "Rfx3 TPM vs Dnahc9 TPM") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RankZ_Rfx3vDnahc9.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

## Correlation between Rfx3 exp vs GFR, male and female separeate##############################################################
Gene_allele <- as.data.frame(Gene_allele)
pheno$Rfx3_allele <- Gene_allele$ENSMUSG00000040929

pheno_F <- pheno[pheno$Sex == "F",]
pheno_F <- pheno_F[complete.cases(pheno_F$C2_log),]
pheno_M <- pheno[pheno$Sex == "M",]
pheno_M <- pheno_M[complete.cases(pheno_M$C2_log),]

fit <- lm(C2_log ~ Rfx3_ranked, pheno_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno_F, aes(y = C2_log, x = Rfx3_ranked)) +
	geom_point(size = 2, aes(colour = Rfx3_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 6.25, x = 0.65, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 1, by = 0.01)) +
	scale_y_continuous( "log-transformed C2 GFR", breaks = seq(0, 6.5, by = 0.5)) +
	labs( title = "RankZ-transformed Rfx3 TPM vs Log-transformed C2 GFR (Females)") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Ranked_Rfx3vC2GFR_Females.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

fit <- lm(C2_log ~ Rfx3_ranked, pheno_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
ggplot <- ggplot(pheno_M, aes(y = C2_log, x = Rfx3_ranked)) +
	geom_point(size = 2, aes(colour = Rfx3_allele)) +
	geom_smooth(method = lm) +
	annotate( "text" , y = 6.25, x = 0.65, label = eq, fontface = "bold", size = 5) + 
	scale_x_continuous( "RankZ-transformed Rfx3 TPM", breaks = seq(0, 1, by = 0.01)) +
	scale_y_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.5)) +
	labs( title = "RanksZ-transformed Rfx3 TPM vs Log-transformed C2 GFR (Males)") +
	theme(plot.title = element_text(hjust = 0.5))  

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Ranked_Rfx3vC2GFR_Males.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

## Rfx3 TPM of founder alleles at Rfx3 locus vs GFR of founder alleles at Rfx3 locus (so basically combining the two plots from this morning

fit <- lm(pheno$Rfx3_tpm ~ pheno$Sex + best.genoprobs.192[,, "UNC30148273"], na.action = na.exclude)
#Female averages by strain at Snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1]  #add female 
#ecfit <- exp(cfit) #exp fixes values back to natural numbers
cfit_F <- as.data.frame(cfit)
#Male averages by strain at snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1] + cfit[2] #add female and male = male (female = female, male = male + female)
#ecfit <- exp(cfit) #exp fixes balues back to natural numbers
cfit_M <- as.data.frame(cfit)
data_coef <- cfit_F
data_coef$M <- cfit_M$cfit
colnames(data_coef) <- c("F","M")
names(data_coef) <- c("F","M")
data_coef <- data_coef[c(3:10),]
rownames(data_coef) <- c("A","B","C","D","E","F","G","H")
data_coef$founder <- rownames(data_coef)
data_Rfx3 <- data_coef

fit <- lm(pheno$C2_log ~ pheno$Sex + best.genoprobs.192[,, "UNC30148273"], na.action = na.exclude)
#Female averages by strain at Snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1]  #add female 
#ecfit <- exp(cfit) #exp fixes values back to natural numbers
cfit_F <- as.data.frame(cfit)
#Male averages by strain at snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1] + cfit[2] #add female and male = male (female = female, male = male + female)
#ecfit <- exp(cfit) #exp fixes balues back to natural numbers
cfit_M <- as.data.frame(cfit)
data_coef <- cfit_F
data_coef$M <- cfit_M$cfit
colnames(data_coef) <- c("F","M")
names(data_coef) <- c("F","M")
data_coef <- data_coef[c(3:10),]
rownames(data_coef) <- c("A","B","C","D","E","F","G","H")
data_coef$founder <- rownames(data_coef)
data_GFR <- data_coef

names(data_Rfx3) <- c("Rfx3_F", "Rfx3_M", "Founder")
names(data_GFR) <- c("GFR_F", "GFR_M", "Founder")

data <- data_Rfx3
data$GFR_F <- data_GFR$GFR_F
data$GFR_M <- data_GFR$GFR_M

##cor by sex
ggplot <- ggplot(data, aes(x = Rfx3_F, y = GFR_F, colour = Founder)) + 
	geom_point(size = 2) +
	scale_x_continuous( "Rfx3 TPM") +
	scale_y_continuous( "log-transformed C2 GFR") +
	labs( title = "Rfx3 TPM vs Log-transformed C2 GFR by Founder(Females)") +
	theme(plot.title = element_text(hjust = 0.5)) 

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3vGFR_Females_by_founder.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

ggplot <- ggplot(data, aes(x = Rfx3_M, y = GFR_M, colour = Founder)) + 
	geom_point(size = 2) +
	scale_x_continuous( "Rfx3 TPM") +
	scale_y_continuous( "log-transformed C2 GFR") +
	labs( title = "Rfx3 TPM vs Log-transformed C2 GFR by Founder(Male)") +
	theme(plot.title = element_text(hjust = 0.5)) 

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Rfx3vGFR_Males_by_founder.pdf", width = 10.0, height = 7.5)
print(ggplot)
dev.off()

## Create narrowing Rfx3 coefficeint plot for grant focusing on LOD peak ###############################################################################
load("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata")

#subset by by position to create a smaller qtl object focusing on target region
qtl <- qtl.GFR.log.C2.192
#subset lod score by postion
lodA <- qtl$lod$A
lodA <- lodA[lodA$chr == 19,]
lodA <- lodA[lodA$pos > 22.5,]
lodA <- lodA[lodA$pos < 35,]
#get corresponding coef scores
coefA <- qtl$coef$A
coefA <- as.data.frame(coefA)
coefA <- coefA[rownames(lodA),]
#reconstruct qtl object for plotting 
qtl$lod$A <- lodA
qtl$coef$A <- coefA

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.GFR.chr19.focused.RFX3.pdf", width = 7, height = 7)
coefplot(qtl, chr = 19, main = "Chr 19 log-transformed C2 GFR Allele Effect Plot")
dev.off()

