library(DOQTL)
library(ggplot2)
library(knitr)
library(reshape2)
library(grid)


setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.ACR10WK.192.Rdata")

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
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

interval = bayesint(qtl.log.ACR10WK.192, chr = 2) #view snp marker
fit <- lm(pheno$ACR10WK_log ~ pheno$Sex + best.genoprobs.192[,,"JAX00098823"], na.action = na.exclude) # rs27466023  2:113483135 

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
	labs(title = "Col4a5xDO JAX00098823/rs27466023 Log-transformed ACR10WK_log by Founder", x = "DO Founders", y = "Log-Transformed ACR10WK") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_JAX00098823_ACR10wkLOG_by_Founders.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

interval = bayesint(qtl.log.Alb10WK.192, chr = 2) #view snp marker
fit <- lm(pheno$Alb10WK_log ~ pheno$Sex + best.genoprobs.192[,,"JAX00098823"], na.action = na.exclude) # rs27466023  2:113483135 

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
	labs(title = "Col4a5xDO JAX00098823/rs27466023 Log-transformed Albumin10WK by Founder", x = "DO Founders", y = "Log-transformed Alb10WK") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Fmn1_JAX00098823_Albumin10wkLOG_by_Founders.pdf", width = 10.0, height = 7.5)	
print(ggplot)
dev.off()

#######




