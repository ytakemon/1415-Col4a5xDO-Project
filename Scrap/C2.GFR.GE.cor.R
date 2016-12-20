#Correlation between log C2 GFR and gene expression
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load sample
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA


GE <- read.delim("/hpcdata/ytakemon/Col4a5xDO/civet_run/160629-134808_1415-1117_cat_R1.fastq.gz/gbrs.quantified.multiway.genes.expected_read_counts", sep = "\t", header = T)
rownames(GE) <- GE[,1]
ACR <- pheno["X1415.1117", "ACR6WK"] #2397.31

temp <- matrix(0, nrow = nrow(GE), ncol = 1)
rownames(temp) <- rownames(GE)
for (i in 1:nrow(temp){
	temp[i,1] <- cor(load()), as.vector(pheno[, "ACR6WK"]))
}




temp[i,x] <- cor(as.vector(Muga[i,,]), as.vector(RNA[x,,]))