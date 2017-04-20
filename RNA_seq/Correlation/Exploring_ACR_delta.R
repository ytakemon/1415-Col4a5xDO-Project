#Exploring delta of ACR
library(reshape2)
library(ggplot2)
library(grid)


# set to working directory on cadillac
setwd("/hpcdata/ytakemon/Col4a5xDO")

#Calculate correlation between phenotype and gene expression via tpm
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
RNA_seq <- as.data.frame(RNA_seq)
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
#clean up pheno and add log of ACR
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
                  "ACR6WK_log", "ACR10WK_log", "ACR15WK_log",
                  "ACR6WK", "ACR10WK", "ACR15WK")]
# Subset out males
pheno_M <- pheno[pheno$Sex == "M",] #95x12

#get delta of ACR15WK to ACR6WK
pheno$delta_ACR15_6 <- pheno$ACR15WK - pheno$ACR6WK
pheno$delta_ACR10_6 <- pheno$ACR10WK - pheno$ACR6WK

pheno_M$delta_ACR15_6 <- pheno_M$ACR15WK - pheno_M$ACR6WK
pheno_M$delta_ACR10_6 <- pheno_M$ACR10WK - pheno_M$ACR6WK

ggdata <- pheno_M[,c("MouseID", "delta_ACR15_6", "delta_ACR10_6")]
ggdata <- melt(ggdata)
ggplot <- ggplot(ggdata, aes(value, ..count.., fill = variable, colour = variable)) +
  geom_density( alpha = 0.1) +
  scale_x_continuous("x axis") +
  scale_y_continuous("y axis") +
  labs( title = "Males ACR delta distribution")+
  coord_flip()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR_deltas_Males.pdf", width = 6, height = 6)
print(ggplot)
dev.off()

ggdata <- pheno[,c("MouseID", "delta_ACR15_6", "delta_ACR10_6")]
ggdata <- melt(ggdata)
ggplot <- ggplot(ggdata, aes(value, ..count.., fill = variable, colour = variable)) +
  geom_density( alpha = 0.1) +
  scale_x_continuous("x axis") +
  scale_y_continuous("y axis") +
  labs( title = "Males ACR delta distribution")+
  coord_flip()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR_deltas_bothsex.pdf", width = 6, height = 6)
print(ggplot)
dev.off()
