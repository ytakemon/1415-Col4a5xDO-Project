#Plotting QTL of the delta between ACR at 15 and 6WKs.
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
pheno$delta_ACR15_6 <- pheno$ACR15WK - pheno$ACR6WK
pheno[pheno ==  -Inf] = NA
options(na.action = "na.pass")
pheno <- pheno[,c("MouseID", "Sex", "ACR6WK", "ACR15WK", "delta_ACR15_6")]

#sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"
sex.covar <- as.matrix(sex.covar)

#GWAS
gwas_delta_ACR15_6 <- scanone.assoc(
	pheno = pheno,
	pheno.col = "delta_ACR15_6",
	probs = best.genoprobs.192,
	K = K_GS,
	addcovar = sex.covar,
	markers = GM_snps,
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz",
	ncl = 3) #number of cores to use
save(gwas_delta_ACR15_6, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.delta_ACR15_6.Rdata")

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS_delta_ACR15_6.pdf", width = 11, height = 7)
plot(gwas_delta_ACR15_6, main = "Col4a5xDO delta ACR15_6wk GWAS map")
dev.off()

##Coefficient plots using QTL Rdata
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl_delta_ACR15_6.192.Rdata")







#ACR
ACR_qtl <- qtl.log.ACR6WK.192
ACR_qtl$coef$A[abs(ACR_qtl$coef$A) > 2 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.ACR6.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(ACR_qtl , chr = 2, main = "Chr 2 ACR6WK_log")
dev.off()
