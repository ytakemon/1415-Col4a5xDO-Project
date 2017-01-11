library(DOQTL)
library(knitr)
library(pcaMethods)
library(RColorBrewer)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#Calculate correlation between phenotype and gene expression via tpm
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
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
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Creat6WK_log","Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log")]
pheno_M <- pheno[pheno$Sex == "M",] #95x9

#GFR
G_pheno <- pheno_M[complete.cases(pheno_M$C2_log),]
G_RNA_seq <- RNA_seq[rownames(G_pheno),]
RNA_GFR_cor <- array(0, c(length(colnames(G_RNA_seq)),1), dimnames = list (colnames(G_RNA_seq), "GFR_corelation"))

for (i in 1 : length(colnames(G_RNA_seq))){
	temp <- cor(G_pheno$C2_log, G_RNA_seq[,i])
	RNA_GFR_cor[i,1] <- temp
}
RNA_GFR_cor <- RNA_GFR_cor[complete.cases(RNA_GFR_cor),] #removes NAs that occured to 0 tpm counts

#Alb6
A6_pheno <- pheno_M[complete.cases(pheno_M$Alb6WK_log),]
A6.RNA_seq <- RNA_seq[rownames(A6_pheno),]
RNA_A6_cor <- array(0, c(length(colnames(A6.RNA_seq)),1), dimnames = list (colnames(A6.RNA_seq), "A6_corelation"))
for (i in 1 : length(colnames(A6.RNA_seq))){
	temp <- cor(A6_pheno$Alb6WK_log, A6.RNA_seq[,i])
	RNA_A6_cor[i,1] <- temp
}
RNA_A6_cor <- RNA_A6_cor[complete.cases(RNA_A6_cor),] #removes NAs that occured to 0 tpm counts


################# Reference



#GFR
G.pheno <- pheno[complete.cases(pheno$C2_log),]
G.RNA_seq <- RNA_seq[rownames(G.pheno),]
RNA_GFR_cor <- array(0, c(length(colnames(G.RNA_seq)),1), dimnames = list (colnames(G.RNA_seq), "GFR_corelation"))

for (i in 1 : length(colnames(G.RNA_seq))){
	temp <- cor(G.pheno$C2_log, G.RNA_seq[,i])
	RNA_GFR_cor[i,1] <- temp
}

RNA_GFR_cor <- RNA_GFR_cor[complete.cases(RNA_GFR_cor),] #removes NAs that occured to 0 tpm counts

#Alb6
A6.pheno <- pheno[complete.cases(pheno$Alb6WK_log),]
A6.RNA_seq <- RNA_seq[rownames(A6.pheno),]
RNA_A6_cor <- array(0, c(length(colnames(A6.RNA_seq)),1), dimnames = list (colnames(A6.RNA_seq), "A6_corelation"))
for (i in 1 : length(colnames(A6.RNA_seq))){
	temp <- cor(A6.pheno$Alb6WK_log, A6.RNA_seq[,i])
	RNA_A6_cor[i,1] <- temp
}
RNA_A6_cor <- RNA_A6_cor[complete.cases(RNA_A6_cor),] #removes NAs that occured to 0 tpm counts
#ACR6
ACR6.pheno <- pheno[complete.cases(pheno$ACR6WK_log),]
ACR6.RNA_seq <- RNA_seq[rownames(ACR6.pheno),]
RNA_ACR6_cor <- array(0, c(length(colnames(ACR6.RNA_seq)),1), dimnames = list (colnames(ACR6.RNA_seq), "ACR6_corelation"))
for (i in 1 : length(colnames(ACR6.RNA_seq))){
	temp <- cor(ACR6.pheno$ACR6WK_log, ACR6.RNA_seq[,i])
	RNA_ACR6_cor[i,1] <- temp
}
RNA_ACR6_cor <- RNA_ACR6_cor[complete.cases(RNA_ACR6_cor),] #removes NAs that occured to 0 tpm counts

#Alb10
A10.pheno <- pheno[complete.cases(pheno$Alb10WK_log),]
A10.RNA_seq <- RNA_seq[rownames(A10.pheno),]
RNA_A10_cor <- array(0, c(length(colnames(A10.RNA_seq)),1), dimnames = list (colnames(A10.RNA_seq), "A10_corelation"))
for (i in 1 : length(colnames(A10.RNA_seq))){
	temp <- cor(A10.pheno$Alb10WK_log, A10.RNA_seq[,i])
	RNA_A10_cor[i,1] <- temp
}
RNA_A10_cor <- RNA_A10_cor[complete.cases(RNA_A10_cor),] #removes NAs that occured to 0 tpm counts
#ACR10
ACR10.pheno <- pheno[complete.cases(pheno$ACR10WK_log),]
ACR10.RNA_seq <- RNA_seq[rownames(ACR10.pheno),]
RNA_ACR10_cor <- array(0, c(length(colnames(ACR10.RNA_seq)),1), dimnames = list (colnames(ACR10.RNA_seq), "ACR10_corelation"))
for (i in 1 : length(colnames(ACR10.RNA_seq))){
	temp <- cor(ACR10.pheno$ACR10WK_log, ACR10.RNA_seq[,i])
	RNA_ACR10_cor[i,1] <- temp
}
RNA_ACR10_cor <- RNA_ACR10_cor[complete.cases(RNA_ACR10_cor),] #removes NAs that occured to 0 tpm counts

#Alb15
A15.pheno <- pheno[complete.cases(pheno$Alb15WK_log),]
A15.RNA_seq <- RNA_seq[rownames(A15.pheno),]
RNA_A15_cor <- array(0, c(length(colnames(A15.RNA_seq)),1), dimnames = list (colnames(A15.RNA_seq), "A15_corelation"))
for (i in 1 : length(colnames(A15.RNA_seq))){
	temp <- cor(A15.pheno$Alb15WK_log, A15.RNA_seq[,i])
	RNA_A15_cor[i,1] <- temp
}
RNA_A15_cor <- RNA_A15_cor[complete.cases(RNA_A15_cor),] #removes NAs that occured to 0 tpm counts
#ACR15
ACR15.pheno <- pheno[complete.cases(pheno$ACR15WK_log),]
ACR15.RNA_seq <- RNA_seq[rownames(ACR15.pheno),]
RNA_ACR15_cor <- array(0, c(length(colnames(ACR15.RNA_seq)),1), dimnames = list (colnames(ACR15.RNA_seq), "ACR15_corelation"))
for (i in 1 : length(colnames(ACR15.RNA_seq))){
	temp <- cor(ACR15.pheno$ACR15WK_log, ACR15.RNA_seq[,i])
	RNA_ACR15_cor[i,1] <- temp
}
RNA_ACR15_cor <- RNA_ACR15_cor[complete.cases(RNA_ACR15_cor),] #removes NAs that occured to 0 tpm counts

max(RNA_GFR_cor)
max(RNA_A6_cor)
max(RNA_A10_cor)
max(RNA_A15_cor)
max(RNA_ACR6_cor)
max(RNA_ACR10_cor)
max(RNA_ACR15_cor)
