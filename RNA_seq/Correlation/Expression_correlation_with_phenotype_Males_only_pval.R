########################################################################################################################
## Expression correlation with phenotype (Males only) with p value
## Author: Yuka Takemon
## Date created: 01/11/17

## OBJECTIVE
##
## sessionInfo()
##
## List of data saved from this script (time savers for reanalysis)
##
########################################################################################################################

# Load libraries
library(DOQTL)

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
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Creat6WK_log",
                  "Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log",
                  "ACR6WK_log", "ACR10WK_log", "ACR15WK_log")]
# Subset out males
pheno_M <- pheno[pheno$Sex == "M",] #95x12

# Phenotpe data has samples where pheotype measurements were not obtained due to QC or other technical errors.
# In order to successfully run cor(), we must remove sampels with NA in phenotype data
# Once that is done, an empty array is created to place the results from cor(gene exp, phenotyp) for each gene.

#GFR
G_pheno <- pheno_M[complete.cases(pheno_M$C2_log),] #remove NAs
G_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(G_pheno),]
RNA_GFR_cor <- array(0, c(length(colnames(G_RNA_seq)),3),
              dimnames = list (colnames(G_RNA_seq), c("GFR_correlation", "df", "pval")))

for (i in 1 : length(colnames(G_RNA_seq))){
	temp <- cor.test(G_pheno$C2_log, G_RNA_seq[,i])
	RNA_GFR_cor[i, "GFR_correlation"] <- temp$estimate[[1]]
  RNA_GFR_cor[i, "df"] <- temp$parameter[[1]]
  RNA_GFR_cor[i, "pval"] <- temp$p.value[[1]]
}
#removes NAs that occured to 0 tpm counts
RNA_GFR_cor <- as.data.frame(RNA_GFR_cor[complete.cases(RNA_GFR_cor),])
RNA_GFR_cor_sig <- RNA_GFR_cor[RNA_GFR_cor$pval < 0.05,]

#ACR6
ACR6_pheno <- pheno_M[complete.cases(pheno_M$ACR6WK_log),]
ACR6_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(ACR6_pheno),]
RNA_ACR6_cor <- array(0, c(length(colnames(ACR6_RNA_seq)),3),
              dimnames = list (colnames(ACR6_RNA_seq), c("ACR6_correlation", "df", "pval")))
for (i in 1 : length(colnames(ACR6_RNA_seq))){
	temp <- cor.test(ACR6_pheno$ACR6WK_log, ACR6_RNA_seq[,i])
	RNA_ACR6_cor[i, "ACR6_correlation"] <- temp$estimate[[1]]
	RNA_ACR6_cor[i, "df"] <- temp$parameter[[1]]
  RNA_ACR6_cor[i, "pval"] <- temp$p.value[[1]]
}
#removes NAs that occured to 0 tpm counts
RNA_ACR6_cor <- as.data.frame(RNA_ACR6_cor[complete.cases(RNA_ACR6_cor),])
RNA_ACR6_cor_sig <- RNA_ACR6_cor[RNA_ACR6_cor$pval < 0.05,]

#ACR10
ACR10_pheno <- pheno_M[complete.cases(pheno_M$ACR10WK_log),]
ACR10_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(ACR10_pheno),]
RNA_ACR10_cor <- array(0, c(length(colnames(ACR10_RNA_seq)),3),
                dimnames = list (colnames(ACR10_RNA_seq), c("ACR10_correlation", "df", "pval")))
for (i in 1 : length(colnames(ACR10_RNA_seq))){
  temp <- cor.test(ACR10_pheno$ACR10WK_log, ACR10_RNA_seq[,i])
  RNA_ACR10_cor[i, "ACR10_correlation"] <- temp$estimate[[1]]
	RNA_ACR10_cor[i, "df"] <- temp$parameter[[1]]
  RNA_ACR10_cor[i, "pval"] <- temp$p.value[[1]]
}
#removes NAs that occured to 0 tpm counts
RNA_ACR10_cor <- as.data.frame(RNA_ACR10_cor[complete.cases(RNA_ACR10_cor),])
RNA_ACR10_cor_sig <- RNA_ACR10_cor[RNA_ACR10_cor$pval < 0.05,]

#ACR15
ACR15_pheno <- pheno_M[complete.cases(pheno_M$ACR15WK_log),]
ACR15_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(ACR15_pheno),]
RNA_ACR15_cor <- array(0, c(length(colnames(ACR15_RNA_seq)),3),
                dimnames = list (colnames(ACR15_RNA_seq), c("ACR15_correlation", "df", "pval")))
for (i in 1 : length(colnames(ACR15_RNA_seq))){
	temp <- cor.test(ACR15_pheno$ACR15WK_log, ACR15_RNA_seq[,i])
  RNA_ACR15_cor[i, "ACR15_correlation"] <- temp$estimate[[1]]
	RNA_ACR15_cor[i, "df"] <- temp$parameter[[1]]
  RNA_ACR15_cor[i, "pval"] <- temp$p.value[[1]]
}
#removes NAs that occured to 0 tpm counts
RNA_ACR15_cor <- as.data.frame(RNA_ACR15_cor[complete.cases(RNA_ACR15_cor),])
RNA_ACR15_cor_sig <- RNA_ACR15_cor[RNA_ACR15_cor$pval < 0.05,]

#To quickly check max and min correlation results
max(RNA_GFR_cor)
max(RNA_ACR6_cor)
max(RNA_ACR10_cor)
max(RNA_ACR15_cor)

min(RNA_GFR_cor)
min(RNA_ACR6_cor)
min(RNA_ACR10_cor)
min(RNA_ACR15_cor)


##Combine biomart data with RNA correlation data
#	Append new columns with corresponding data.
#	Sort data first by ensembl gene id before subsetting and appending to the list, otherwise there might be errors with subsetting

#load Ensemble ID form biomart (preivously saved on cadillac)
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")

# there are duplicates in the Ensembl_mouseID, but chromosome name and start and end postons
# are the same so i am only extracting the unique samples
Ensembl_mouseID <- Ensembl_mouseID[!duplicated(Ensembl_mouseID$ensembl_gene_id),]
rownames(Ensembl_mouseID) <- make.names(Ensembl_mouseID[,1])
Ensembl_mouseID <- Ensembl_mouseID[ order( rownames(Ensembl_mouseID)),]

#subset Ensembl_mouseID to match cor data.frame
RNA_GFR_cor$geneID <- rownames(RNA_GFR_cor)
RNA_GFR_cor <- RNA_GFR_cor[ order( rownames(RNA_GFR_cor)),]
Ensembl_mouseID_G <- Ensembl_mouseID[rownames(RNA_GFR_cor),]

RNA_ACR6_cor$geneID <- rownames(RNA_ACR6_cor)
RNA_ACR6_cor <- RNA_ACR6_cor[ order( rownames(RNA_ACR6_cor)),]
Ensembl_mouseID_ACR6 <- Ensembl_mouseID[rownames(RNA_ACR6_cor),]

RNA_ACR10_cor$geneID <- rownames(RNA_ACR10_cor)
RNA_ACR10_cor <- RNA_ACR10_cor[ order( rownames(RNA_ACR10_cor)),]
Ensembl_mouseID_ACR10 <- Ensembl_mouseID[rownames(RNA_ACR10_cor),]

RNA_ACR15_cor$geneID <- rownames(RNA_ACR15_cor)
RNA_ACR15_cor <- RNA_ACR15_cor[ order( rownames(RNA_ACR15_cor)),]
Ensembl_mouseID_ACR15 <- Ensembl_mouseID[rownames(RNA_ACR15_cor),]

###NOTES: NAs are introduced in rownames of because those ENS genes were removed from the current ensembl database, due to lack of
###   additional evidence
#	Examples of removed ensemble ids (Look up on ensembl website for most up to date info)
#	"ENSMUSG00000103861"
#	"ENSMUSG00000103888"
#	"ENSMUSG00000104003"

# Append ensembl information to correlation data for human readable format.
RNA_GFR_cor$hgnc_symbol <- Ensembl_mouseID_G$hgnc_symbol
RNA_GFR_cor$mgi_symbol <- Ensembl_mouseID_G$mgi_symbol
RNA_GFR_cor$chromosome <- Ensembl_mouseID_G$chromosome_name
RNA_GFR_cor$start <- Ensembl_mouseID_G$start_position
RNA_GFR_cor$end <- Ensembl_mouseID_G$end_position

RNA_ACR6_cor$hgnc_symbol <- Ensembl_mouseID_ACR6$hgnc_symbol
RNA_ACR6_cor$mgi_symbol <- Ensembl_mouseID_ACR6$mgi_symbol
RNA_ACR6_cor$chromosome <- Ensembl_mouseID_ACR6$chromosome_name
RNA_ACR6_cor$start <- Ensembl_mouseID_ACR6$start_position
RNA_ACR6_cor$end <- Ensembl_mouseID_ACR6$end_position

RNA_ACR10_cor$hgnc_symbol <- Ensembl_mouseID_ACR10$hgnc_symbol
RNA_ACR10_cor$mgi_symbol <- Ensembl_mouseID_ACR10$mgi_symbol
RNA_ACR10_cor$chromosome <- Ensembl_mouseID_ACR10$chromosome_name
RNA_ACR10_cor$start <- Ensembl_mouseID_ACR10$start_position
RNA_ACR10_cor$end <- Ensembl_mouseID_ACR10$end_position

RNA_ACR15_cor$hgnc_symbol <- Ensembl_mouseID_ACR15$hgnc_symbol
RNA_ACR15_cor$mgi_symbol <- Ensembl_mouseID_ACR15$mgi_symbol
RNA_ACR15_cor$chromosome <- Ensembl_mouseID_ACR15$chromosome_name
RNA_ACR15_cor$start <- Ensembl_mouseID_ACR15$start_position
RNA_ACR15_cor$end <- Ensembl_mouseID_ACR15$end_position

#Order correlation in decreasing order
RNA_GFR_cor <- RNA_GFR_cor[ order( RNA_GFR_cor$GFR_correlation, decreasing = TRUE),]
RNA_ACR6_cor <- RNA_ACR6_cor[ order( RNA_ACR6_cor$ACR6_correlation, decreasing = TRUE),]
RNA_ACR10_cor <- RNA_ACR10_cor[ order( RNA_ACR10_cor$ACR10_correlation, decreasing = TRUE),]
RNA_ACR15_cor <- RNA_ACR15_cor[ order( RNA_ACR15_cor$ACR15_correlation, decreasing = TRUE),]

# Extract only necessary columns to show.
RNA_GFR_cor <- RNA_GFR_cor[,c("geneID", "GFR_correlation", "df", "pval", "mgi_symbol", "chromosome", "start", "end")]
RNA_GFR_cor_sig <- RNA_GFR_cor[RNA_GFR_cor$pval < 0.05,]
RNA_ACR6_cor <- RNA_ACR6_cor[,c("geneID", "ACR6_correlation", "df", "pval", "mgi_symbol", "chromosome", "start", "end")]
RNA_ACR6_cor_sig <- RNA_ACR6_cor[RNA_ACR6_cor$pval < 0.05,]
RNA_ACR10_cor <- RNA_ACR10_cor[,c("geneID", "ACR10_correlation", "df", "pval", "mgi_symbol", "chromosome", "start", "end")]
RNA_ACR10_cor_sig <- RNA_ACR10_cor[RNA_ACR10_cor$pval < 0.05,]
RNA_ACR15_cor <- RNA_ACR15_cor[,c("geneID", "ACR15_correlation", "df", "pval", "mgi_symbol", "chromosome", "start", "end")]
RNA_ACR15_cor_sig <- RNA_ACR15_cor[RNA_ACR15_cor$pval < 0.05,]

# Export as tab sep txt files, and remove rownames because that displaces column names by one.
write.table(RNA_GFR_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_GFR_M_cor.txt", sep = "\t", row.names = FALSE)
write.table(RNA_A6_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_Alb6_M_cor.txt", sep = "\t", row.names = FALSE)
write.table(RNA_A10_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_Alb10_M_cor.txt", sep = "\t", row.names = FALSE)
write.table(RNA_A15_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_Alb15_M_cor.txt", sep = "\t", row.names = FALSE)
write.table(RNA_ACR6_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_ACR6_M_cor.txt", sep = "\t", row.names = FALSE)
write.table(RNA_ACR10_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_ACR10_M_cor.txt", sep = "\t", row.names = FALSE)
write.table(RNA_ACR15_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_ACR15_M_cor.txt", sep = "\t", row.names = FALSE)

##################
#MAKE SCATTER PLOT OF ACR AT 15WK AND GFR FOR MALES
#load data
gfr <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_GFR_M_cor.txt",
					sep = "\t", header = TRUE)

acr15 <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_ACR15_M_cor.txt",
					sep = "\t", header = TRUE)
#subset data
rownames(acr15) <- acr15$geneID
rownames(gfr) <- gfr$geneID

#> dim(acr15)
#[1] 35332     8
#> dim(gfr)
#[1] 34565     8

#subset out genes that are presetn in both acr15 and gfr
sub_acr <- acr15[acr15$geneID %in% rownames(gfr),]
sub_acr <- sub_acr[ order(rownames(sub_acr)),]

sub_gfr <- gfr[gfr$geneID %in% rownames(sub_acr),]
sub_gfr <- sub_gfr[ order(rownames(sub_gfr)),]
#> dim(sub_acr)
#[1] 34396     8
#> dim(sub_gfr)
#[1] 34396     8

data <- sub_acr
data$RNA_GFR_cor <- sub_gfr$RNA_GFR_cor
data <- data[,c("geneID", "RNA_ACR15_cor", "RNA_GFR_cor", "mgi_symbol", "chromosome", "start", "end")]

write.table(data, "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_pheno_data/males/RNA_GFR_ACR_M_merge_cor.txt", sep = "\t")

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_GFR_ACR_M_merge.pdf", width = 10.0, height = 7.5)
plot(data$RNA_GFR_cor, data$RNA_ACR15_cor)
dev.off()
