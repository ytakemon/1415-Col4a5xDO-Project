#1415 Col4a5xDO qtl 
#Yuka Takemon
#11/02/16

## Description:
#This script randomizes sample within the sexes and set set as "male.perm" and "female.perm" for each permutation
#The set randomization is then assigned to Hf and Ha QTL analysis. 
#The Hi is calculated by Hi = Hf - Ha, but only for the LOD scores
#All other values of Hi is "fake". Hi object is a shell which retains DOQTL formatting to hold unique LOD scores for compatibility

## Notes:
#Scanone (qtl) permutations have not yet implemented permutations for interactive covariates, so it must be done manually
#this script only extracts from Autosomes
#input args to specify how many perms to do
#I think the limit to qsub per user is around 500, so will need to split this in to 1:500, 501:1000 submissions

##############################
#R script for looping
library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")
args <- commandArgs(trailingOnly = TRUE)

#load files
#genoprobs
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
#GM_grid
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
#kinship
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
#phenotype
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean up phenotype data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#add log of GFR
pheno[pheno < 0 ] = NA
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2)
options(na.action = 'na.pass') #leave in NAs

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#need to randomize samples by scrambling all the females and males separately
#subset sexes
males <- which(sex.covar[,1] == 1)
females <- which(sex.covar[,1] == 0)

#randomize samples using sample()
males.perm <- sample( rownames( pheno)[males]) 
females.perm <- sample( rownames(pheno)[females])

#repace rownames of pheno data with set randomized order of samples
rownames(pheno)[males] <- males.perm 
rownames(pheno)[females] <- females.perm
#repace rownames of sex.covar data with set randomized order of samples
rownames(sex.covar)[males] <- males.perm
rownames(sex.covar)[females] <- females.perm

#set as matrix to feed into qlt
sex.covar <- as.matrix(sex.covar)

#run qtl 
Hf_qtl <- scanone(pheno = pheno, pheno.col = "C2_log", probs =  best.genoprobs.192, K = K_GS, addcovar = sex.covar, intcovar = sex.covar, snps = GM_snps)
Ha_qtl <- scanone(pheno = pheno, pheno.col = "C2_log", probs =  best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)

#Create shell for Hi and calculate LODi score
Hi_qtl <- Hf_qtl #creates shell
Hi_qtl$lod$A$lod <- Hi_qtl$lod$A$lod - Ha_qtl$lod$A$lod #replaces LOD score with LODf - LODa

#find max lod score for A
max.lod <- max(Hi_qtl$lod$A$lod)
write( max.lod, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.Hi.perm/lod.perm.",args[1],".txt"))

#End of Loop


############################
#This portion is in linux Cadillac
cd /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/GFR.Hi.perm/
grep -h . *.txt > 1000.perms.txt
#check number of lines/ permutation in file
wc -l 1000.perms.txt

###########################
#back in R interactive processing with: QTL_sex_intcovar_GFR.R




