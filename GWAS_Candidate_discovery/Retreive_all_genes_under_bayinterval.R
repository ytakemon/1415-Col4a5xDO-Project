#get gene list of every interval above the 0.6 green line in the QTL map
library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load sample
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)
pheno$delta_ACR15_6 <- (pheno$ACR15WK - pheno$ACR6WK)

pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Alb10WK_log <- log(pheno$Alb10WK)
pheno$Alb15WK_log <- log(pheno$Alb15WK)

pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno$Creat10WK_log <- log(pheno$Creat10WK)
pheno$Creat15WK_log <- log(pheno$Creat15WK)

pheno$C2_log <- log(pheno$C2)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#create covar with sex and creatinine
#All 6wk
temp <- sex.covar
temp$creat6wk <- pheno$Creat6WK_log
covar.sex.creat6wk <- temp
#All 10wk
temp <- sex.covar
temp$creat10wk <- pheno$Creat10WK_log
covar.sex.creat10wk <- temp
#All 15wk
temp <- sex.covar
temp$creat15wk <- pheno$Creat15WK_log
covar.sex.creat15wk <- temp

sex.covar <- as.matrix(sex.covar)
covar.sex.creat6wk <- as.matrix(covar.sex.creat6wk)
covar.sex.creat10wk <- as.matrix(covar.sex.creat10wk)
covar.sex.creat15wk <- as.matrix(covar.sex.creat15wk)
#kinship mapping
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")

#################################################################################
# Gather list of genes using UCSC list
UCSC_list <- read.delim("/hpcdata/ytakemon/genome_index/Mus_musculus/UCSC/mm10/GenomeStudio/Mus_musculus/UCSC-mm10/refGene.txt", header = FALSE)
UCSC_list <- UCSC_list[,c(13,2,3,5,6)]
names(UCSC_list) <- c("Gene_name", "UCSC_id", "Chr", "Start", "End")

#Starting with C2 GFR gene list becuase its less overwhelming
GFR_gene_list <- data.frame(Gene_name=character(),
                 UCSC_id=character(),
                 Chr=integer(),
                 Start=integer(),
                 End=integer(),
                 stringsAsFactors=FALSE)
#Chr2
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 20000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 42000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 56000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 65000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

#Chr4
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr4",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 0,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 10000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

#Chr7
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr7",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 120000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 126000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

#Chr9
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr9",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 30000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 36000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

#Chr10
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr10",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 88000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 96000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

#Chr15
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr15",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 50000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 55000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr15",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 86000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 89000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

#Chr19
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr19",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 25000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 29000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
GFR_gene_list <- rbind(GFR_gene_list, sub_UCSC_list)

rownames(GFR_gene_list) <- NULL
write.table(GFR_gene_list, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/GFR_gene_list.txt", sep = "\t")

# Albumin6wk_wCreatinine_covar_list
Alb6WK_gene_list <- data.frame(Gene_name=character(),
                 UCSC_id=character(),
                 Chr=integer(),
                 Start=integer(),
                 End=integer(),
                 stringsAsFactors=FALSE)

#Chr2
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 111000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 114000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 172000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 175000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr4
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr4",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 120000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 128000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr5
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr5",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 105000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 116000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr7
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr7",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 115000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 122000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr9
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr9",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 119000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 122000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr11
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr11",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 4000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 10000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr11",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 71000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 76000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr11",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 87000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 90000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr12
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr12",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 75000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 82000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr12",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 108000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 110000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

#Chr16
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr12",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 90000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 93000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb6WK_gene_list <- rbind(Alb6WK_gene_list, sub_UCSC_list)

rownames(Alb6WK_gene_list) <- NULL
write.table(Alb6WK_gene_list, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/Alb6WK_log_gene_list.txt", sep = "\t")

# Albumin10wk_wCreatinine_covar_list
Alb10WK_gene_list <- data.frame(Gene_name=character(),
                 UCSC_id=character(),
                 Chr=integer(),
                 Start=integer(),
                 End=integer(),
                 stringsAsFactors=FALSE)

#Chr2
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 108000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 115000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb10WK_gene_list <- rbind(Alb10WK_gene_list, sub_UCSC_list)

#Chr4
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr4",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 62000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 65000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb10WK_gene_list <- rbind(Alb10WK_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr4",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 120000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 125000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb10WK_gene_list <- rbind(Alb10WK_gene_list, sub_UCSC_list)

#Chr11
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr11",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 0,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 10000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb10WK_gene_list <- rbind(Alb10WK_gene_list, sub_UCSC_list)

#Chr13
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr13",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 100000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 105000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb10WK_gene_list <- rbind(Alb10WK_gene_list, sub_UCSC_list)

rownames(Alb10WK_gene_list) <- NULL
write.table(Alb10WK_gene_list, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/Alb10WK_log_gene_list.txt", sep = "\t")

# Albumin15wk_wCreatinine_covar_list
Alb15WK_gene_list <- data.frame(Gene_name=character(),
                 UCSC_id=character(),
                 Chr=integer(),
                 Start=integer(),
                 End=integer(),
                 stringsAsFactors=FALSE)

#Chr15
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr15",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 86000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 89000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
Alb15WK_gene_list <- rbind(Alb15WK_gene_list, sub_UCSC_list)

rownames(Alb15WK_gene_list) <- NULL
write.table(Alb15WK_gene_list, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/Alb15WK_log_gene_list.txt", sep = "\t")

# delta_ACR15_6_list
delta_ACR15_6_gene_list <- data.frame(Gene_name=character(),
                 UCSC_id=character(),
                 Chr=integer(),
                 Start=integer(),
                 End=integer(),
                 stringsAsFactors=FALSE)

#Chr1
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr1",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 36000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 42000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

#Chr2
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 35000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 43000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 80000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 92000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 110000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 125000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

#Chr5
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr5",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 49000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 54000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr5",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 65000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 74000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

#Ch13
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "ch13",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 39000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 42000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

#Chr14
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr14",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 80000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 90000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

#Chr17
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr17",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 25000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 30000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr17",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 66000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 73000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

#Chr19
sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr19",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 20000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 25000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$Gene_name),]
delta_ACR15_6_gene_list <- rbind(delta_ACR15_6_gene_list, sub_UCSC_list)

rownames(delta_ACR15_6_gene_list) <- NULL
write.table(delta_ACR15_6_gene_list, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/delta_ACR15_6_gene_list.txt", sep = "\t")

# load gene list files again
GFR_gene_list <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/GFR_gene_list.txt",
                            sep = "\t", header = TRUE)
Alb6WK_gene_list <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/Alb6WK_log_gene_list.txt",
                            sep = "\t", header = TRUE)
Alb10WK_gene_list <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/Alb10WK_log_gene_list.txt",
                            sep = "\t", header = TRUE)
Alb15WK_gene_list <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/Alb15WK_log_gene_list.txt",
                            sep = "\t", header = TRUE)
delta_ACR15_6_gene_list <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/delta_ACR15_6_gene_list.txt",
                            sep = "\t", header = TRUE)
# Get gene name only
GFR_genes <- as.character(GFR_gene_list$Gene_name)
Alb6WK_genes <- as.character(Alb6WK_gene_list$Gene_name)
Alb10WK_genes <- as.character(Alb10WK_gene_list$Gene_name)
Alb15WK_genes <- as.character(Alb15WK_gene_list$Gene_name)
delta_ACR15_6_genes <- as.character(delta_ACR15_6_gene_list$Gene_name)

#Get intersect
GFR_Alb6WK_intersect <- intersect(GFR_genes, Alb6WK_genes)
GFR_Alb10WK_intersect <- intersect(GFR_genes, Alb10WK_genes)
GFR_Alb15WK_intersect <- intersect(GFR_genes, Alb15WK_genes)
GFR_delta_ACR15_6_intersect <- intersect(GFR_genes, delta_ACR15_6_genes)

Alb6WK_Alb10WK_intersect <- intersect(Alb6WK_genes, Alb10WK_genes)
Alb6WK_Alb15WK_intersect <- intersect(Alb6WK_genes, Alb15WK_genes)
Alb6WK_delta_ACR15_6_intersect <- intersect(Alb6WK_genes, delta_ACR15_6_genes)

Alb10WK_Alb15WK_intersect <- intersect(Alb10WK_genes, Alb15WK_genes)
Alb10WK_delta_ACR15_6_intersect <- intersect(Alb10WK_genes, delta_ACR15_6_genes)

Alb15WK_delta_ACR15_6_intersect <- intersect(Alb15WK_genes, delta_ACR15_6_genes)

# Apend all into one list
list <- list(GFR_Alb6WK_intersect = GFR_Alb6WK_intersect, GFR_Alb10WK_intersect = GFR_Alb10WK_intersect, GFR_Alb15WK_intersect = GFR_Alb15WK_intersect, GFR_delta_ACR15_6_intersect =GFR_delta_ACR15_6_intersect,
            Alb6WK_Alb10WK_intersect = Alb6WK_Alb10WK_intersect, Alb6WK_Alb15WK_intersect = Alb6WK_Alb15WK_intersect, Alb6WK_delta_ACR15_6_intersect = Alb6WK_delta_ACR15_6_intersect,
            Alb10WK_Alb15WK_intersect = Alb10WK_Alb15WK_intersect, Alb10WK_delta_ACR15_6_intersect = Alb10WK_delta_ACR15_6_intersect,
            Alb15WK_delta_ACR15_6_intersect = Alb15WK_delta_ACR15_6_intersect)

lapply(list, function(x) write.table( data.frame(x), './GBRS_reconstruction/reconstruct/best.compiled.genoprob/compiled_gene_list/All_intersect_genes.txt'  , append= TRUE, sep='\t' ))
