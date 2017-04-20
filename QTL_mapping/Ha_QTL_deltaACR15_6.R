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

#run QTL
#qtl <- scanone( pheno = pheno, pheno.col = "delta_ACR15_6",
#                probs = best.genoprobs.192, K = K_GS,
#                addcovar = sex.covar, snps = GM_snps)
#save(qtl, file = "GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl_delta_ACR15_6.192.Rdata")

#permutations are run via Ha_QTL_deltaACR15_6_perms.R
load("GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl_delta_ACR15_6.192.Rdata")



pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Ha_QTL_deltaACR15_6_perms.pdf", width = 10.0, height = 7.5)
plot(qtl, main = "Col4a5xDO 192 genoprob delta_ACR15_6")
dev.off()
