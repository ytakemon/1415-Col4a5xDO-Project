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
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA

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

#need to randomize samples by scrambling all the females and males separately
#subset sexes
males <- which(sex.covar[,1] == 1)
females <- which(sex.covar[,1] == 0)

#randomize samples using sample()
males.perm <- sample( rownames( pheno)[males]) 
females.perm <- sample( rownames(pheno)[females])

#replace rownames of pheno data with set randomized order of samples
rownames(pheno)[males] <- males.perm 
rownames(pheno)[females] <- females.perm
#replace rownames of sex.covar data with set randomized order of samples
rownames(sex.covar)[males] <- males.perm
rownames(sex.covar)[females] <- females.perm
#replace rownames of covar.sex.creat
rownames(covar.sex.creat6wk)[males] <- males.perm
rownames(covar.sex.creat6wk)[females] <- females.perm

rownames(covar.sex.creat10wk)[males] <- males.perm
rownames(covar.sex.creat10wk)[females] <- females.perm

rownames(covar.sex.creat15wk)[males] <- males.perm
rownames(covar.sex.creat15wk)[females] <- females.perm

#set as matrix to feed into qlt
sex.covar <- as.matrix(sex.covar)
covar.sex.creat6wk <- as.matrix(covar.sex.creat6wk)
covar.sex.creat10wk <- as.matrix(covar.sex.creat10wk)
covar.sex.creat15wk <- as.matrix(covar.sex.creat15wk)

#run qtl with randomized samples
Hf_6 <- scanone(pheno = pheno, pheno.col = "Alb6WK_log", probs =  best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat6wk, intcovar = sex.covar, snps = GM_snps)
Hf_10 <- scanone(pheno = pheno, pheno.col = "Alb10WK_log", probs =  best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat10wk, intcovar = sex.covar, snps = GM_snps)
Hf_15 <- scanone(pheno = pheno, pheno.col = "Alb15WK_log", probs =  best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat15wk, intcovar = sex.covar, snps = GM_snps)

Ha_6 <- scanone(pheno = pheno, pheno.col = "Alb6WK_log", probs =  best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat6wk, snps = GM_snps)
Ha_10 <- scanone(pheno = pheno, pheno.col = "Alb10WK_log", probs =  best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat10wk, snps = GM_snps)
Ha_15 <- scanone(pheno = pheno, pheno.col = "Alb15WK_log", probs =  best.genoprobs.192, K = K_GS, addcovar = covar.sex.creat15wk, snps = GM_snps)

#Create shell for Hi and calculate LODi score
Hi_6 <- Hf_6
Hi_10 <- Hf_10
Hi_15 <- Hf_15

#Calculate and replace LOD score by LODf - LODa
Hi_6$lod$A$lod <- Hi_6$lod$A$lod - Ha_6$lod$A$lod
Hi_10$lod$A$lod <- Hi_10$lod$A$lod - Ha_10$lod$A$lod
Hi_15$lod$A$lod <- Hi_15$lod$A$lod - Ha_15$lod$A$lod

max.lod_6 <- max(Hi_6$lod$A$lod)
max.lod_10 <- max(Hi_10$lod$A$lod)
max.lod_15 <- max(Hi_15$lod$A$lod)

write( max.lod_6, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.Hi.perm/lod.perm.",args[1],".txt"))
write( max.lod_10, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.Hi.perm/lod.perm.",args[1],".txt"))
write( max.lod_15, file = paste0("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.Hi.perm/lod.perm.",args[1],".txt"))

#End of Loop

############################
#This portion is in linux Cadillac
cd /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb6.Hi.perm/
grep -h . *.txt > 1000.perms.txt
#check number of lines/ permutation in file
wc -l 1000.perms.txt

cd /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb10.Hi.perm/
grep -h . *.txt > 1000.perms.txt
#check number of lines/ permutation in file
wc -l 1000.perms.txt

cd /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl.perm/Alb15.Hi.perm/
grep -h . *.txt > 1000.perms.txt
#check number of lines/ permutation in file
wc -l 1000.perms.txt

###########################
#back in R interactive processing with: QTL_sex_intcovar_Alb.R










