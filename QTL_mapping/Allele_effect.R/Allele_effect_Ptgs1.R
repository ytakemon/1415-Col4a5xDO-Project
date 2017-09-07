library(DOQTL)
library(ggplot2)
library(knitr)
library(reshape2)
library(grid)

setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR6WK.192.Rdata")

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs

# Plcd1 is location Chromosome 9: 119,071,527-119,093,502
# Get list of markers in this position then find marker with highest lod score
qtl <- qtl.log.ACR6WK.192
qtl <- qtl$lod$A
qtl <- qtl[qtl$chr == 9,]
qtl <- qtl[qtl$pos > 119.071,]
qtl <- qtl[qtl$pos < 119.094,]
target <- qtl[qtl$lod == max(qtl$lod),] # pos 101.761

#Look up rsid
GM_snps[GM_snps$marker == target$marker,] # rs29691525

#calculate coef at target marker
fit <- lm(pheno$ACR6WK_log ~ pheno$Sex + best.genoprobs.192[,,target$marker], na.action = na.exclude)

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
	labs(title = paste0("Col4a5xDO allele effect at ", target$marker, " position: ", target$pos, " for Log-transformed ACR6WK by founder"),
				x = "DO Founders",
				y = "Log-Transformed ACR6WK values") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

#Female averages by strain at Snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1]  #add female
ecfit <- exp(cfit) #exp fixes values back to natural numbers
ecfit_F <- as.data.frame(ecfit)
#Male averages by strain at snp marker
cfit <- coef(fit)
cfit[10] = 0 #set standard (ie strain compared to) to 0
cfit <- cfit + cfit[1] + cfit[2] #add female and male = male (female = female, male = male + female)
ecfit <- exp(cfit) #exp fixes balues back to natural numbers
ecfit_M <- as.data.frame(ecfit)

data_coef <- ecfit_F
data_coef$M <- ecfit_M$ecfit
colnames(data_coef) <- c("F","M")
names(data_coef) <- c("F","M")
data_coef <- data_coef[c(3:10),]
rownames(data_coef) <- c("A","B","C","D","E","F","G","H")
data_coef$founder <- rownames(data_coef)

eggdata <- melt(data_coef)
colnames(eggdata) <- c("Founder", "Sex", "Value")
eggplot <- ggplot(eggdata, aes( x = Founder, y = Value, fill = Sex)) +
	geom_bar( stat = "identity", position = "dodge") +
	scale_fill_discrete( labels = c("Females","Males")) +
	labs(title = paste0("Col4a5xDO allele effect at Plcd1 (chr", target$chr, " ", target$marker, " position: ", target$pos, ") for ACR6WK by founder"),
				x = "DO Founders",
				y = "ACR6WK values") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

pdf(paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Allele_effect_at_", target$marker,"_ACR6WK_Plcd1.pdf"), width = 10.0, height = 7.5)
print(eggplot)
dev.off()
