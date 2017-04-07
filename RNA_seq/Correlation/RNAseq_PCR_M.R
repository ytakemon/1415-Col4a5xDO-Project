# Load libraries
library(DOQTL)
library(pcaMethods)
library(grid)
library(ggplot2)

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
# Subset out males
pheno_M <- pheno[pheno$Sex == "M",] #95x12
RNA_seq_M <- RNA_seq[rownames(RNA_seq) %in% rownames(pheno_M),]

#PCA analysis row = sample column = variables
pca_data <- pca(object = RNA_seq_M, nPcs = nrow(RNA_seq_M))
pca_score <- as.data.frame(scores(pca_data))



ggplot12 <- ggplot(pca_score, aes(x = PC1, y =PC2))+
			geom_point()
ggplot32 <- ggplot(pca_score, aes(x = PC3, y =PC2))+
			geom_point()
ggplot34 <- ggplot(pca_score, aes(x = PC3, y =PC4))+
			geom_point()
ggplot54 <- ggplot(pca_score, aes(x = PC5, y =PC4))+
			geom_point()


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PCA_M_pc12345.pdf", width = 6, height = 6)
pushViewport(viewport( layout = grid.layout(2,2)))
print(ggplot12, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot32, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ggplot34, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(ggplot54, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()



data_pca <- pca( object = t(data), nPcs = ncol(data))
data_pca_score <- test <- as.data.frame(scores(data_pca))
data_pca_score$geno <- "WT"
data_pca_score$geno[4:6] <- "KO"

pdf("/hpcdata/ytakemon/Far2KO/Plots/pca.pdf")
ggplot(data_pca_score, aes(x= PC1, y = PC2, colour = geno)) +
geom_point()
dev.off()
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a5_fig2_ACR_GFR_dist_YUKA.pdf", width = 12, height = 6)
pushViewport(viewport(layout = grid.layout(1, 2)))
print(ggplot_F, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot_M, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()