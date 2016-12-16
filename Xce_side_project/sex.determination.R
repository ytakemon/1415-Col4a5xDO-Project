## Sex determination of 4 samples that appears to be mis- phenotyped reavearled by PCA plots
#	Samples that phenotyped M, but appeared F :
#X1415.0914
#X1415.0243
#	Samples that phenotyped F, but appeared M :
#X1415.1091
#X1415.1016

#Sex speficif gene expression to check sex
#Jarid1c for X specific ENSMUSG00000025332
#Jarid1d for Y specific ENSMUSG00000056673

library(DOQTL)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_tpm.Rdata")
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
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Alb10WK_log","Alb15WK_log","ACR6WK_log", "ACR10WK_log", "ACR15WK_log")]

#check to see Col4 genes are present
geneID <- colnames(RNA_seq)
geneID[geneID == "ENSMUSG00000025332"] #Jarid1c aka Kdm5c (found)
geneID[geneID == "ENSMUSG00000056673"] #Jarid1d aka Kdm5d (found)
#Lookin good

#Extract genes into 1 data frame
Col_names <- c("Kdm5c_X_tpm", "Kdm5d_Y_tpm")
Sex_RNA_tpm <- array(0, c(192, 2), dimnames = list(rownames(pheno), Col_names))
for (i in 1:192){

	Sex_RNA_tpm[i, "Kdm5c_X_tpm"] <- RNA_seq[i, "ENSMUSG00000025332"]
	Sex_RNA_tpm[i, "Kdm5d_Y_tpm"] <- RNA_seq[i, "ENSMUSG00000056673"]

}
Sex_RNA_tpm <- as.data.frame(Sex_RNA_tpm)
Sex_RNA_tpm$Sex_pheno <- pheno$Sex

#The samples that were phenotyped as males did have male specific expression, however their expression seems significantly lower than the average
#The samples that were phenotyped as females did have female specific expression, however their expression seems significantly lower than the average

#subset samples and do T-test
Female_all <- Sex_RNA_tpm[Sex_RNA_tpm$Sex_pheno == "F",] 
Female_all$issue <- "Female"
Female_all["X1415.1016", "issue"] <- "Problem"
Female_all["X1415.1091", "issue"] <- "Problem"

Male_all <- Sex_RNA_tpm[Sex_RNA_tpm$Sex_pheno == "M",] 
Male_all$issue <- "Male"
Male_all["X1415.0243", "issue"] <- "Problem"
Male_all["X1415.0914", "issue"] <- "Problem"

Female <- ggplot( Female_all, aes(x  = issue, y = Kdm5c_X_tpm)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =15)) +
	xlab( "Female") +
	ylab( "Kdm5c_X_tpm") +
	ggtitle("Female expression of Kdm5c_X_tpm")+
	theme(plot.title = element_text(hjust = 0.5)) 

Male_Y <- ggplot( Male_all, aes(x  = issue, y = Kdm5d_Y_tpm)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =15)) +
	xlab( "Male") +
	ylab( "Kdm5d_Y_tpm") +
	ggtitle("Male expression of Kdm5d_Y_tpm")+
	theme(plot.title = element_text(hjust = 0.5)) 

Male_X <- ggplot( Male_all, aes(x  = issue, y = Kdm5c_X_tpm)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =15)) +
	xlab( "Male") +
	ylab( "Kdm5c_X_tpm") +
	ggtitle("Male expression of Kdm5c_X_tpm")+
	theme(plot.title = element_text(hjust = 0.5)) 

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Sex_expression.png", width = 1000, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 2)))
print(Female, vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1))
print(Male_X, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(Male_Y, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()








