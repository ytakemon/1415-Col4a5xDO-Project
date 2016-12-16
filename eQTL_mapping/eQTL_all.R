## eQTL analysis of all transcripts from RNA-seq data
## Yuka Takemon
## created: 11/30/16
## Modified last: 12/13/16

#	This script will create a scanone qtl objct as well as a QTL map for all genes from RNA-seq give an numeric input into arugument.
#	The intenet is to run on cadillac given a loop job submission 1 - 46,517 genes. However due to Cadillac's user restrictions, we can
#	only submit 1000 jobs on queue at one time.

#	resources requested to complet each task: nodes=1:ppn=1,walltime=00:20:00

#load library
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#Input data
args <- commandArgs(trailingOnly = TRUE) ## input for looped script in bash
args <- as.numeric(args)
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#outputdir
qtl_outdir <- "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/"
plot_outdir <- "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Complete_eQTL_plots/"

RNA_seqZ <- as.data.frame(RNA_seqZ)

Ensembl_mouseID <- Ensembl_mouseID[!duplicated(Ensembl_mouseID$ensembl_gene_id),]
rownames(Ensembl_mouseID) <- make.names(Ensembl_mouseID[,1]) 
Ensembl_mouseID <- Ensembl_mouseID[ order( rownames(Ensembl_mouseID)),]
Ensembl_mouseID <- Ensembl_mouseID[rownames(t(RNA_seqZ)),]

#create sex covariate
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#Get list of genes from RNA_seq
RNA_list <- colnames(RNA_seqZ) # there are  46517 genes

#qtl
temp <- scanone( pheno = RNA_seqZ, pheno.col = RNA_list[args[1]], probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
object_name <- paste(RNA_list[args[1]], Ensembl_mouseID[args[1], "mgi_symbol"], "eQTL", sep = "." )
assign(object_name, temp)
save( list = object_name, file = paste( qtl_outdir, RNA_list[args[1]], ".", 
	Ensembl_mouseID[args[1],"mgi_symbol"], ".eQTL.rankZ.tpm.Rdata", sep = ""))

#plot qtl
pdf( paste(plot_outdir, RNA_list[args[1]],".", Ensembl_mouseID[args[1],"mgi_symbol"],".eQTL.perm1000.pdf", sep = ""), 
	width = 10.0, height = 7.5)
plot(temp, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = paste("eQTL plot of ",  RNA_list[args[1]],
	" ", Ensembl_mouseID[args[1],"mgi_symbol"] ," w/ 1000 perms", sep = ""))
dev.off()


#################################################################################################################################
## There were some interruptions when running the above script on Cadillac, which resulted in a few unknown genes not being
## run. We have a total of 46,405 genes out of the 46,517 genes run. We need to first identify the 112 genes that did not get run 
## and figure out what # they are on the correpsonding list to run.

## Helper functions:
#	Subset string from the right, x = character string, n = number of characters from the right
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
#	Subset string from the left, x = character string, n = number of characters from the left
substrLeft <- function(x, n){
  substr(x, 1, n)
}
## Do first before starting this section:
#	Run everyting above besides scanone and plotting. 

#	Get list of file that were completed previously
All_eqtl_files <- list.files("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata", full.names = T)

#	Extract first 18 characters of the QLT file string, standard ENSEMBL ID has 18 characters (7 alpha and 11 numeric, eg. ENSMUSG00000040929)
#	First, remove common directory path, then extract the first 18 charactes 
All_eqtl_files <- sub("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/","",All_eqtl_files)
All_eqtl_files <- substrLeft(All_eqtl_files, 18)
#	subset out the 112 ENSEMBL IDs that are not in ALL_eqtl_files object.
missing <- setdiff(RNA_list, All_eqtl_files) #element in RNA_list not in All_eqtl_files
save(missing, file ="./GBRS_reconstruction/reconstruct/best.compiled.genoprob/missing.112.genes.Rdata")
##################################################################################################################################
#altered loop script (All_eQTL.R) to run missing 112 genes on cadillac:
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#load data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/missing.112.genes.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#Define argument input 
args <- commandArgs(trailingOnly = TRUE) ## input for looped script in bash
args <- as.numeric(args)
args <- missing[args]

#outputdir
qtl_outdir <- "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/"
plot_outdir <- "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Complete_eQTL_plots/"

RNA_seqZ <- as.data.frame(RNA_seqZ)

Ensembl_mouseID <- Ensembl_mouseID[!duplicated(Ensembl_mouseID$ensembl_gene_id),]
rownames(Ensembl_mouseID) <- make.names(Ensembl_mouseID[,1]) 
Ensembl_mouseID <- Ensembl_mouseID[ order( rownames(Ensembl_mouseID)),]
Ensembl_mouseID <- Ensembl_mouseID[rownames(t(RNA_seqZ)),]

#create sex covariate
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

temp <- scanone( pheno = RNA_seqZ, pheno.col = args, probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
object_name <- paste(args, Ensembl_mouseID[args[1], "mgi_symbol"], "eQTL", sep = "." )
assign(object_name, temp)
save( list = object_name, file = paste( qtl_outdir, args, ".", 
	Ensembl_mouseID[args[1],"mgi_symbol"], ".eQTL.rankZ.tpm.Rdata", sep = ""))

#plot qtl
pdf( paste(plot_outdir, args,".", Ensembl_mouseID[args[1],"mgi_symbol"],".eQTL.perm1000.pdf", sep = ""), 
	width = 10.0, height = 7.5)
plot(temp, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = paste("eQTL plot of ",  args,
	" ", Ensembl_mouseID[args[1],"mgi_symbol"] ," w/ 1000 perms", sep = ""))
dev.off()







