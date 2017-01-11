########################################################################################################################
## Expression correlation with phenotype 
## Date created: 12/11/16 

## OBJECTIVE
## We want to see if there are any significant correlation between gene expression and phenotype that may suggest dysregulation
## of certian pathways and mechanisms. As this is the first analysis using RNA-seq gene expression data, we will need to 
## run a QC with PC analysis to make sure there aren't any obvious underlying factors that effects the global patterning of 
## gene expression. 

## sessionInfo()
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## locale:
## [1] C
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets
## [8] methods   base
## other attached packages:
## [1] RColorBrewer_1.1-2   knitr_1.11           DOQTL_1.0.0
## [4] AnnotationDbi_1.28.2 GenomeInfoDb_1.2.5   IRanges_2.0.1
## [7] S4Vectors_0.4.0      Biobase_2.26.0       BiocGenerics_0.12.1
## [10] RSQLite_1.0.0        DBI_0.3.1
## loaded via a namespace (and not attached):
##  [1] Biostrings_2.34.1      GenomicRanges_1.18.4   MUGAExampleData_1.0.0
##  [4] QTLRel_0.2-14          RCurl_1.95-4.7         RUnit_0.4.30
##  [7] Rsamtools_1.18.3       XML_3.98-1.3           XVector_0.6.0
## [10] annotate_1.44.0        annotationTools_1.40.0 biomaRt_2.22.0
## [13] bitops_1.0-6           corpcor_1.6.8          gdata_2.17.0
## [16] gtools_3.5.0           hwriter_1.3.2          mclust_5.1
## [19] org.Hs.eg.db_3.0.0     org.Mm.eg.db_3.0.0     tools_3.1.1
## [22] xtable_1.8-0           zlibbioc_1.12.0

## List of data saved from this script (time savers for reanalysis)
## Rdata:
## save(RNA_seq, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
## save(mouseID, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")
## Tables:
## write.table(RNA_GFR_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_GFR_cor.txt", sep = "\t", row.names = TRUE)
## write.table(RNA_A6_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_Alb6_cor.txt", sep = "\t", row.names = TRUE)
## write.table(RNA_A10_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_Alb10_cor.txt", sep = "\t", row.names = TRUE)
## write.table(RNA_A15_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_Alb15_cor.txt", sep = "\t", row.names = TRUE)
## write.table(RNA_ACR6_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_ACR6_cor.txt", sep = "\t", row.names = TRUE)
## write.table(RNA_ACR10_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_ACR10_cor.txt", sep = "\t", row.names = TRUE)
## write.table(RNA_ACR15_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_ACR15_cor.txt", sep = "\t", row.names = TRUE)
## Plots:
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC_percvar.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC_plot1.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC_plot1_outlier.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC1234_plot.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GFR_cor.png", width = 1500, height = 1000, res = 
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GFR_rank.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A6_cor.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A6_rank.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A10_cor.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A10_rank.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A15_cor.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A15_rank.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6_cor.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6_rank.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10_cor.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10_rank.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR15_cor.png", width = 1500, height = 1000, res = 100)
## png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR15_rank.png", width = 1500, height = 1000, res = 100)

########################################################################################################################

#load necessary data
library(DOQTL)
library(knitr)

# set to working directory
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load essential data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
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

#Consolidate and create RNA_seq data table 
civet_dir <- list.files("./civet_run/", full.names = T) # Extract full paths to each alignment
civet_sample_list <- read.delim("./Sample_list/RNAseq_1415_sample_list.txt", header = FALSE) # Get sample names formatted as civet 1415-0000
genoprob_names <- rownames(best.genoprobs.192) # extract sample names we have genoprobs data for.
Complete_set_names <- as.data.frame(civet_dir) # reformat as data.frame
Complete_set_names$civet_names <- civet_sample_list # Give sample names in submitted format (1415-XXXX) for each sample path
Complete_set_names$genoprob_names <- genoprob_names # Give more consistant names for each sampels (1415_XXXX)
colnames(Complete_set_names) <- names(Complete_set_names) # not sure why colnames and names are out of sync, but force it in sync here. 

#create temp file from tpm data from civet output so I can exract its rownames and column names to be used i lin 65 when consolidating results
temp <- read.delim(file = paste(Complete_set_names[1, "civet_dir"],"/", "gbrs.quantified.diploid.genes.tpm", sep =""), header = TRUE, sep ="\t")
rownames(temp) <- temp$locus
temp$locus <- NULL
col_name <- colnames(temp)
locus_names <- rownames(temp)

#all tpm files are named the samed housed under a directory named after the sample. 
#need to extract files and name them after the sample for distinction
#compiling into single array like genoprobs dataset
RNA_seq <- array(0, c(192, 46517), dimnames = list (Complete_set_names$genoprob_names, locus_names))
for (i in 1:192){
	temp <- read.delim(file = paste(Complete_set_names[i, "civet_dir"],"/", "gbrs.quantified.diploid.genes.tpm", sep =""), header = TRUE, sep ="\t")
	RNA_seq[i,] <- temp$total
}
save(RNA_seq, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
#The file is now saved, so I dont have to redo this again due to extensive processing time.


#RNA_counts <- array(0, c(192, 46517), dimnames = list (Complete_set_names$genoprob_names, locus_names))
#for (i in 1:192){
#	temp <- read.delim(file = paste(Complete_set_names[i, "civet_dir"],"/", "gbrs.quantified.diploid.genes.expected_read_counts", sep =""), header = TRUE, sep ="\t")
#	RNA_counts[i,] <- temp$total
#}
#save(RNA_counts, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_expected_counts.Rdata")

### PCA plots
tRNA_seq <- t(RNA_seq)
expr <- apply(tRNA_seq, 2, rankZ)
pcexpr <- pca(object = t(expr), nPcs = ncol(expr))

perc.var <- sDev(pcexpr)^2 / sum(sDev(pcexpr)^2) * 100
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC_percvar.png", width = 1500, height = 1000, res = 100)
plot(perc.var, main = "Percent variaion of all PC1-192")
dev.off()

cols <- as.numeric(factor(pheno$Sex))
pch <- 16 

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC_plot1.png", width = 1500, height = 1000, res = 100)
plot(scores(pcexpr), col = cols, pch = pch, main = "PC1 vs PC2 (male =red, felame = black)", )
textxy(X = scores(pcexpr)[,1], Y = scores(pcexpr)[,2], labs = rownames(scores(pcexpr)))
dev.off()

#subset to identify outlier 
tmp <- as.data.frame(scores(pcexpr))
tmp <- tmp[tmp$PC1 > 0 ,]
tmp <- tmp[tmp$PC2 < -3 ,]

pheno_sub <- pheno[rownames(tmp),]
cols <- as.numeric(factor(pheno_sub$Sex))

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC_plot1_outlier.png", width = 1500, height = 1000, res = 100)
plot(tmp[,1:2], col = cols, pch = pch. main = "PC1 vs PC2 outlier close up")
textxy(X = tmp[,1], Y = tmp[,2], labs = rownames(tmp))
dev.off()
##outliers identified:
#	X1415.0914 currntly labled M
#	X1415.0243 currntly labled M
#	X1415.1091 currntly labled F
#	X1415.1016 currntly labled F
#Holly is doing a PCR to check from whole tail.

cols <- as.numeric(factor(pheno$Sex))
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PC1234_plot.png", width = 1500, height = 1000, res = 100)
layout(matrix(1:4,2,2))
plot(scores(pcexpr)[,1:2], pch = pch, col = cols)
plot(scores(pcexpr)[,5:4], pch = pch, col = cols)
plot(scores(pcexpr)[,3:2], pch = pch, col = cols)
plot(scores(pcexpr)[,3:4], pch = pch, col = cols)
dev.off()

#######
#extra things for options
#xce <- as.numeric(factor(Xce_allele$Xce_winner))
#pch <- 14+xce
#pch[pch == 16] <- 8
#total_counts <- rowSums(RNA_counts)
#total_cor <- cor(total_counts, scores(pcexpr))
#######

#Calculate correlation between phenotype and gene expression via tpm
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_expected_counts.Rdata")

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

##BIOMART to gather ensembl id to append to RNA_correlation data
#	This ony needs to be done onceto get the latest version of ID and correspond gene symbols that are more meaningful to us.
#	Will also contain data such as chromosome and start and end site of the gene.

# use biomaRt to append a column with gene name
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt) #cannot have DOQTL package loaded for this. So either restart R environment or detach DOQTL and all of its depenednencies 
#current build:mmusculus_gene_ensembl Mus musculus genes (GRCm38.p4) GRCm38.p4
#mart <- useEnsembl("ENSEMBL_MART_ENSEMBL")
#listDatasets(mart) #dataset to use  mmusculus_gene_ensembl,
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", verbose = TRUE)
#head(listAttributes(ensembl)) #get list of attributes to input to getBM
mouseID <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "mgi_symbol", "chromosome_name", "start_position", "end_position"), mart = ensembl)
save(mouseID, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")

##Combine biomart data with RNA correlation data
#	append new columns with corresponding data.
#	probably data needs to get sorted first by ensembl gene id before subsetting and appending to the list
RNA_GFR_cor <- as.data.frame(RNA_GFR_cor)
RNA_A6_cor <- as.data.frame(RNA_A6_cor)
RNA_A10_cor <- as.data.frame(RNA_A10_cor)
RNA_A15_cor <- as.data.frame(RNA_A15_cor)
RNA_ACR6_cor <- as.data.frame(RNA_ACR6_cor)
RNA_ACR10_cor <- as.data.frame(RNA_ACR10_cor)
RNA_ACR15_cor <- as.data.frame(RNA_ACR15_cor)
mouseID <- as.data.frame(mouseID)

#there are duplicates in the mouseID, but chromosome name and start and end postons are the same so i am only extracting the unique sampels
mouseID <- mouseID[!duplicated(mouseID$ensembl_gene_id),]
rownames(mouseID) <- make.names(mouseID[,1]) 
mouseID <- mouseID[ order( rownames(mouseID)),]

#subset mouseID to match cor data
#> dim(mouseID)
#[1] 49100     6
#> dim(RNA_GFR_cor)
#[1] 35766     1
RNA_GFR_cor$geneID <- rownames(RNA_GFR_cor)
RNA_GFR_cor <- RNA_GFR_cor[ order( rownames(RNA_GFR_cor)),]
mouseID_G <- mouseID[rownames(RNA_GFR_cor),]
#> dim(mouseID_G)
#[1] 35766     6

#> dim(RNA_A6_cor)
#[1] 36539     1
RNA_A6_cor$geneID <- rownames(RNA_A6_cor)
RNA_A6_cor <- RNA_A6_cor[ order( rownames(RNA_A6_cor)),]
mouseID_A6 <- mouseID[rownames(RNA_A6_cor),]
#> dim(mouseID_A6)
#[1] 36539     6

#> dim(RNA_A10_cor)
#[1] 36480     1
RNA_A10_cor$geneID <- rownames(RNA_A10_cor)
RNA_A10_cor <- RNA_A10_cor[ order( rownames(RNA_A10_cor)),]
mouseID_A10 <- mouseID[rownames(RNA_A10_cor),]
#> dim(mouseID_A10)
#[1] 36480     6

#> dim(RNA_A15_cor)
#[1] 36554     1
RNA_A15_cor$geneID <- rownames(RNA_A15_cor)
RNA_A15_cor <- RNA_A15_cor[ order( rownames(RNA_A15_cor)),]
mouseID_A15 <- mouseID[rownames(RNA_A15_cor),]
#> dim(mouseID_A15)
#[1] 36554     6

#> dim(RNA_ACR6_cor)
#[1] 36333     1
RNA_ACR6_cor$geneID <- rownames(RNA_ACR6_cor)
RNA_ACR6_cor <- RNA_ACR6_cor[ order( rownames(RNA_ACR6_cor)),]
mouseID_ACR6 <- mouseID[rownames(RNA_ACR6_cor),]
#> dim(mouseID_ACR6)
#[1] 36333     6

#> dim(RNA_ACR10_cor)
#[1] 36321     1
RNA_ACR10_cor$geneID <- rownames(RNA_ACR10_cor)
RNA_ACR10_cor <- RNA_ACR10_cor[ order( rownames(RNA_ACR10_cor)),]
mouseID_ACR10 <- mouseID[rownames(RNA_ACR10_cor),]
#> dim(mouseID_ACR10)
#[1] 36321     6

#> dim(RNA_ACR15_cor)
#[1] 36436     1
RNA_ACR15_cor$geneID <- rownames(RNA_ACR15_cor)
RNA_ACR15_cor <- RNA_ACR15_cor[ order( rownames(RNA_ACR15_cor)),]
mouseID_ACR15 <- mouseID[rownames(RNA_ACR15_cor),]
#> dim(mouseID_ACR15)
#[1] 36436     6

###NOTES: NAs are introduced in rowname because those ENS genes were removed from the data based and is no longer in the current database
#	Examples:
#	"ENSMUSG00000103861"
#	"ENSMUSG00000103888"
#	"ENSMUSG00000104003"

RNA_GFR_cor$hgnc_symbol <- mouseID_G$hgnc_symbol
RNA_GFR_cor$mgi_symbol <- mouseID_G$mgi_symbol
RNA_GFR_cor$chromosome <- mouseID_G$chromosome_name
RNA_GFR_cor$start <- mouseID_G$start_position
RNA_GFR_cor$end <- mouseID_G$end_position

RNA_A6_cor$hgnc_symbol <- mouseID_A6$hgnc_symbol
RNA_A6_cor$mgi_symbol <- mouseID_A6$mgi_symbol
RNA_A6_cor$chromosome <- mouseID_A6$chromosome_name
RNA_A6_cor$start <- mouseID_A6$start_position
RNA_A6_cor$end <- mouseID_A6$end_position

RNA_A10_cor$hgnc_symbol <- mouseID_A10$hgnc_symbol
RNA_A10_cor$mgi_symbol <- mouseID_A10$mgi_symbol
RNA_A10_cor$chromosome <- mouseID_A10$chromosome_name
RNA_A10_cor$start <- mouseID_A10$start_position
RNA_A10_cor$end <- mouseID_A10$end_position

RNA_A15_cor$hgnc_symbol <- mouseID_A15$hgnc_symbol
RNA_A15_cor$mgi_symbol <- mouseID_A15$mgi_symbol
RNA_A15_cor$chromosome <- mouseID_A15$chromosome_name
RNA_A15_cor$start <- mouseID_A15$start_position
RNA_A15_cor$end <- mouseID_A15$end_position

RNA_ACR6_cor$hgnc_symbol <- mouseID_ACR6$hgnc_symbol
RNA_ACR6_cor$mgi_symbol <- mouseID_ACR6$mgi_symbol
RNA_ACR6_cor$chromosome <- mouseID_ACR6$chromosome_name
RNA_ACR6_cor$start <- mouseID_ACR6$start_position
RNA_ACR6_cor$end <- mouseID_ACR6$end_position

RNA_ACR10_cor$hgnc_symbol <- mouseID_ACR10$hgnc_symbol
RNA_ACR10_cor$mgi_symbol <- mouseID_ACR10$mgi_symbol
RNA_ACR10_cor$chromosome <- mouseID_ACR10$chromosome_name
RNA_ACR10_cor$start <- mouseID_ACR10$start_position
RNA_ACR10_cor$end <- mouseID_ACR10$end_position

RNA_ACR15_cor$hgnc_symbol <- mouseID_ACR15$hgnc_symbol
RNA_ACR15_cor$mgi_symbol <- mouseID_ACR15$mgi_symbol
RNA_ACR15_cor$chromosome <- mouseID_ACR15$chromosome_name
RNA_ACR15_cor$start <- mouseID_ACR15$start_position
RNA_ACR15_cor$end <- mouseID_ACR15$end_position

RNA_GFR_cor <- RNA_GFR_cor[ order( RNA_GFR_cor$RNA_GFR_cor, decreasing = TRUE),]
RNA_A6_cor <- RNA_A6_cor[ order( RNA_A6_cor$RNA_A6_cor, decreasing = TRUE),]
RNA_A10_cor <- RNA_A10_cor[ order( RNA_A10_cor$RNA_A10_cor, decreasing = TRUE),]
RNA_A15_cor <- RNA_A15_cor[ order( RNA_A15_cor$RNA_A15_cor, decreasing = TRUE),]
RNA_ACR6_cor <- RNA_ACR6_cor[ order( RNA_ACR6_cor$RNA_ACR6_cor, decreasing = TRUE),]
RNA_ACR10_cor <- RNA_ACR10_cor[ order( RNA_ACR10_cor$RNA_ACR10_cor, decreasing = TRUE),]
RNA_ACR15_cor <- RNA_ACR15_cor[ order( RNA_ACR15_cor$RNA_ACR15_cor, decreasing = TRUE),]

#rankZ transform tpm
G_rank <- apply(as.matrix(RNA_GFR_cor$RNA_GFR_cor), 2, rankZ) #need library(DOQTL)
A6_rank <- apply(as.matrix(RNA_A6_cor$RNA_A6_cor), 2, rankZ)
A10_rank <- apply(as.matrix(RNA_A10_cor$RNA_A10_cor), 2, rankZ)
A15_rank <- apply(as.matrix(RNA_A15_cor$RNA_A15_cor), 2, rankZ)
ACR6_rank <- apply(as.matrix(RNA_ACR6_cor$RNA_ACR6_cor), 2, rankZ)
ACR10_rank <- apply(as.matrix(RNA_ACR10_cor$RNA_ACR10_cor), 2, rankZ)
ACR15_rank <- apply(as.matrix(RNA_ACR15_cor$RNA_ACR15_cor), 2, rankZ)

RNA_GFR_cor$rankZ <- G_rank
RNA_A6_cor$rankZ <- A6_rank
RNA_A10_cor$rankZ <- A10_rank
RNA_A15_cor$rankZ <- A15_rank
RNA_ACR6_cor$rankZ <- ACR6_rank
RNA_ACR10_cor$rankZ <- ACR10_rank
RNA_ACR15_cor$rankZ <- ACR15_rank

RNA_GFR_cor <- RNA_GFR_cor[,c("geneID", "RNA_GFR_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]
RNA_A6_cor <- RNA_A6_cor[,c("geneID", "RNA_A6_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]
RNA_A10_cor <- RNA_A10_cor[,c("geneID", "RNA_A10_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]
RNA_A15_cor <- RNA_A15_cor[,c("geneID", "RNA_A15_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]
RNA_ACR6_cor <- RNA_ACR6_cor[,c("geneID", "RNA_ACR6_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]
RNA_ACR10_cor <- RNA_ACR10_cor[,c("geneID", "RNA_ACR10_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]
RNA_ACR15_cor <- RNA_ACR15_cor[,c("geneID", "RNA_ACR15_cor", "rankZ", "hgnc_symbol", "mgi_symbol", "chromosome", "start", "end")]

write.table(RNA_GFR_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_GFR_cor.txt", sep = "\t", row.names = TRUE)
write.table(RNA_A6_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_Alb6_cor.txt", sep = "\t", row.names = TRUE)
write.table(RNA_A10_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_Alb10_cor.txt", sep = "\t", row.names = TRUE)
write.table(RNA_A15_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_Alb15_cor.txt", sep = "\t", row.names = TRUE)
write.table(RNA_ACR6_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_ACR6_cor.txt", sep = "\t", row.names = TRUE)
write.table(RNA_ACR10_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_ACR10_cor.txt", sep = "\t", row.names = TRUE)
write.table(RNA_ACR15_cor, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_ACR15_cor.txt", sep = "\t", row.names = TRUE)
####tests
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GFR_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_GFR_cor$RNA_GFR_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GFR_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_GFR_cor$rankZ)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A6_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_A6_cor$RNA_A6_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A6_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_A6_cor$rankZ)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A10_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_A10_cor$RNA_A10_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A10_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_A10_cor$rankZ)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A15_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_A15_cor$RNA_A15_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/A15_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_A15_cor$rankZ)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_ACR6_cor$RNA_ACR6_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR6_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_ACR6_cor$rankZ)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_ACR10_cor$RNA_ACR10_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR10_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_ACR10_cor$rankZ)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR15_cor.png", width = 1500, height = 1000, res = 100)
plot(RNA_ACR15_cor$RNA_ACR15_cor)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/ACR15_rank.png", width = 1500, height = 1000, res = 100)
plot(RNA_ACR15_cor$rankZ)
dev.off()






















