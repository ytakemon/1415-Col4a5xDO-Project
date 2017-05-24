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
#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#create a subset separating Male and male samples
best.genoprobs.192 <- best.genoprobs.192[ order( rownames( best.genoprobs.192)),,]
pheno <- pheno[ order( rownames(pheno)),]

pheno.F <- subset(pheno, Sex =="F")
pheno.M <- subset(pheno, Sex == "M")

genoprob.F <- best.genoprobs.192[rownames(pheno.F),,]
genoprob.M <- best.genoprobs.192[rownames(pheno.M),,]

sex.covar$MouseID <- rownames(sex.covar)
sex.covar <- sex.covar[ order( sex.covar$MouseID),]

sex.covar.F <- sex.covar[rownames(pheno.F),]
sex.covar.M <- sex.covar[rownames(pheno.M),]

sex.covar$MouseID <- NULL
sex.covar.F$MouseID <- NULL
sex.covar.M$MouseID <- NULL

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

# Alb6wk Gwas Chr2 genes udner Alb10wk bayesian interval.
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno,
	pheno.col ="Alb6WK_log",
	probs = best.genoprobs.192,
	K = K_GS[[chr]],
	addcovar = covar.sex.creat6wk,
	snps = GM_snps,
	chr = 2,
	start = 101,
	end = 114,
	output = "p-value")

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb6.Chr2.10wkinterval.candidates.pdf", width = 10.0, height = 7.5)
assoc.plot(chr2.genes, thr = 5, show.sdps = TRUE)
dev.off()

#Get list of gens between pos 101 - 114  101,000,067 - 113,999,957
UCSC_list <- read.delim("/hpcdata/ytakemon/genome_index/Mus_musculus/UCSC/mm10/GenomeStudio/Mus_musculus/UCSC-mm10/refGene.txt", header = FALSE)
UCSC_list <- UCSC_list[,c(13,2,3,5,6)]
names(UCSC_list) <- c("Gene_name", "UCSC_id", "Chr", "Start", "End")

sub_UCSC_list <- UCSC_list[UCSC_list$Chr == "chr2",]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$Start > 101000000,]
sub_UCSC_list <- sub_UCSC_list[sub_UCSC_list$End < 114000000,]
sub_UCSC_list <- sub_UCSC_list[order(sub_UCSC_list$Start),]
sub_UCSC_list <- sub_UCSC_list[!duplicated(sub_UCSC_list$UCSC_id),]

write.table(sub_UCSC_list, "/hpcdata/ytakemon/Col4a5xDO/Chr2_interval_gene_list.txt", sep = "\t")

#Get list with Ensembl genes names
library(rtracklayer) # to easily parse gtf files
gtf <- readGFF("/hpcdata/ytakemon/genome_index/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf")
names(gtf)[1] <- "Chr"
genes <- gtf[,c("gene_id", "gene_name", "Chr","start", "end")]
sub_Ensembl_list <- genes[genes$Chr == 2,]
sub_Ensembl_list <- sub_Ensembl_list[sub_Ensembl_list$start > 101000000,]
sub_Ensembl_list <- sub_Ensembl_list[sub_Ensembl_list$end < 114000000,]
sub_Ensembl_list <- sub_Ensembl_list[order(sub_Ensembl_list$start),]
sub_Ensembl_list <- sub_Ensembl_list[!duplicated(sub_Ensembl_list$gene_id),]
rownames(sub_Ensembl_list) <- NULL
