#source and laod DOQTL package: this only needs to happen once
#source("http://bioconductor.org/biocLite.R")
#biocLite("DOQTL")
#library(devtools)
#library(knitr)
#install_github("dmgatti/DOQTL", force = TRUE)


library(DOQTL)
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/")

# Load data -------------------------------------------------------------------------------------------

#load geno probs data
load("~/Desktop/Col4a5xDO/R_col4a5xdo/new_col4a5_probs.Rdata")

#create kinship probablility from JAX ftp site containg GigaMuga array snps. 
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#remove I for 9th genotype, as this will mess up the 8 founder IDs
colnames(probs) = sub("I", "", colnames(probs)) 

# Clean up data -----------------------------------------------------------------------------------

# Set the sample IDs for probs to all be 4 digits long, consistent with pheno data
spl = strsplit(rownames(probs), "\\.")
one = sapply(spl, "[", 1)
two = sapply(spl, "[", 2)
wh = which(nchar(two) == 3)
two[wh] = paste0("0", two[wh])
rownames(probs) = paste(one, two, sep = ".")

#6 samples were identified to be B6 samples 75-86, the problem hasn't been identified
#one samples 31 was duplicated with 24
#all 7 were removed since the problem currnetly could not be resolved
probs = probs[-c(24,31,75,80,81,84,85,86),,] #end result 175 samples left

#load in phenotype data containing sex, generation, GFR, and ACR
pheno <- read.table("../Phenotype/1415_master_pheno.txt", header = TRUE, sep= "\t")
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names

pheno[pheno < 0 ] = NA

#07/20/16 202 samples were sent to GeneSeek, 174 samples were successfully genotyped with good reconstruction, total 28 samples needs repeating
# pheno contains 206 samples, 
#dim(pheno)
#dim(probs)

#extract samples that exist in both pheno and probs should be 174 samples
samples = intersect(rownames(pheno), rownames(probs))
pheno = pheno[samples, ]
probs = probs[samples,,]

#run kinship matrix
K_GM = kinship.probs(probs, snps = GM_snps, bychr = TRUE)

##Kinship mapping
image(1:nrow(K_GM[[1]]), 1:ncol(K_GM[[1]]), K_GM[[1]][,ncol(K_GM[[1]]):1], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship of 1415 Col4a5xDO GM", 
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 20 * 0:7, labels = 20 * 7:0, las = 1)

#create covarates
options(na.action = 'na.pass') #leave in NAs
#GFR_pheno = model.matrix(~0+GFR_C2, data = pheno) # 0+ removes intercept column, GFR covar
#GFR_pheno[ GFR_pheno < 0 ] = NA
#GFR_pheno <- na.omit(GFR_pheno)
#GFR_pheno <- data.matrix(GFR_pheno, rownames.force = NA)
#ACR6_pheno = model.matrix(~0+ACR6, data = pheno)
#ACR10_pheno = model.matrix(~0+ACR10, data = pheno)
#ACR15_pheno = model.matrix(~0+ACR15, data = pheno)
Sex_covar = model.matrix(~0+Sex, data = pheno)
colnames(Sex_covar)[2] = "sex"
Sex_covar <- Sex_covar[,"sex"] #removes extra columns now sex, F=0, M=1
Sex_covar <- as.data.frame(Sex_covar)
colnames(Sex_covar)[1] = "sex"

#check distrubution---------------------------------------------------------------------------------------
#ACR
tiff("./Plots/ACR/ACR6WK_distribution.tiff")
hist( pheno$ACR6WK, main = paste("Distribution of 6wk ACR"), xlab = "ACR", ylab = "Frequency")
dev.off()

tiff("./Plots/ACR/ACR10WK_distribution.tiff")
hist( pheno$ACR10WK, main = paste("Distribution of 10wk ACR"), xlab = "ACR", ylab = "Frequency")
dev.off()

tiff("./Plots/ACR/ACR15WK_distribution.tiff")
hist( pheno$ACR15WK, main = paste("Distribution of 15wk ACR"), xlab = "ACR", ylab = "Frequency")
dev.off()

#logACR
pheno$LogACR6WK = log(pheno$ACR6WK)
pheno$LogACR10WK = log(pheno$ACR10WK)
pheno$LogACR15WK = log(pheno$ACR15WK)

tiff("./Plots/logACR/LogACR6WK_distribution.tiff")
hist( pheno$LogACR6WK, main = paste("Distribution of 6wk log[ACR]"), xlab = "log[ACR]", ylab = "Frequency")
dev.off()

tiff("./Plots/logACR/LogACR10WK_distribution.tiff")
hist( pheno$LogACR10WK, main = paste("Distribution of 10wk log[ACR]"), xlab = "log[ACR]", ylab = "Frequency")
dev.off()

tiff("./Plots/logACR/LogACR15WK_distribution.tiff")
hist( pheno$LogACR15WK, main = paste("Distribution of 15wk log[ACR]"), xlab = "log[ACR]", ylab = "Frequency")
dev.off()


#Scanone for ACR qtl at 6wks, 10wks, 15wks -----------------------------------------------------------------------
ACR_6wks_qtl = scanone( pheno = pheno, pheno.col = "ACR6WK", probs = probs, K = K_GM,
                   addcovar = Sex_covar, snps = GM_snps)
ACR_10wks_qtl = scanone( pheno = pheno, pheno.col = "ACR10WK", probs = probs, K = K_GM,
                        addcovar = Sex_covar, snps = GM_snps)
ACR_15wks_qtl = scanone( pheno = pheno, pheno.col = "ACR15WK", probs = probs, K = K_GM,
                        addcovar = Sex_covar, snps = GM_snps)

Log_ACR_6wks_qtl = scanone( pheno = pheno, pheno.col = "LogACR6WK", probs = probs, K = K_GM,
                        addcovar = Sex_covar, snps = GM_snps)
Log_ACR_10wks_qtl = scanone( pheno = pheno, pheno.col = "LogACR10WK", probs = probs, K = K_GM,
                         addcovar = Sex_covar, snps = GM_snps)
Log_ACR_15wks_qtl = scanone( pheno = pheno, pheno.col = "LogACR15WK", probs = probs, K = K_GM,
                         addcovar = Sex_covar, snps = GM_snps)

#plot qtl
tiff("./Plots/ACR/6wk_ACR_QTL.tiff")
plot(ACR_6wks_qtl, main ="Col4a5xDO 6wk ACR QTL")
dev.off()
tiff("./Plots/ACR/10wk_ACR_QTL.tiff")
plot(ACR_10wks_qtl, main ="Col4a5xDO 10wk ACR QTL")
dev.off()
tiff("./Plots/ACR/15wk_ACR_QTL.tiff")
plot(ACR_15wks_qtl, main ="Col4a5xDO 15wk ACR QTL")
dev.off()

tiff("./Plots/LogACR/6wk_LogACR_QTL.tiff")
plot(Log_ACR_6wks_qtl, main ="Col4a5xDO 6wk Log[ACR] QTL")
dev.off()
tiff("./Plots/LogACR/10wk_LogACR_QTL.tiff")
plot(Log_ACR_10wks_qtl, main ="Col4a5xDO 10wk log[ACR] QTL")
dev.off()
tiff("./Plots/LogACR/15wk_LogACR_QTL.tiff")
plot(Log_ACR_15wks_qtl, main ="Col4a5xDO 15wk Log[ACR] QTL")
dev.off()

# Create permuations for significante for plotting ---------------------------------------------------------------

#cerate permuations for significance testing.
perms_6wkACR_1000 <- scanone.perm(pheno = pheno, pheno.col = "ACR6WK", probs = probs,
                                 addcovar = Sex_covar, snps = GM_snps, nperm = 1000)
perms_10wkACR_1000 <- scanone.perm(pheno = pheno, pheno.col = "ACR10WK", probs = probs,
                                 addcovar = Sex_covar, snps = GM_snps, nperm = 1000)
perms_15wkACR_1000 <- scanone.perm(pheno = pheno, pheno.col = "ACR15WK", probs = probs,
                                 addcovar = Sex_covar, snps = GM_snps, nperm = 1000)

perms_Log6wkACR_1000 <- scanone.perm(pheno = pheno, pheno.col = "LogACR6WK", probs = probs,
                                  addcovar = Sex_covar, snps = GM_snps, nperm = 1000)
perms_Log10wkACR_1000 <- scanone.perm(pheno = pheno, pheno.col = "LogACR10WK", probs = probs,
                                   addcovar = Sex_covar, snps = GM_snps, nperm = 1000)
perms_Log15wkACR_1000 <- scanone.perm(pheno = pheno, pheno.col = "LogACR15WK", probs = probs,
                                   addcovar = Sex_covar, snps = GM_snps, nperm = 1000)

#create threshold for alphas
thr_6wkACR_1000 <- get.sig.thr( perms_6wkACR_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr_10wkACR_1000 <- get.sig.thr( perms_10wkACR_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr_15wkACR_1000 <- get.sig.thr( perms_15wkACR_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

thr_Log6wkACR_1000 <- get.sig.thr( perms_Log6wkACR_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr_Log10wkACR_1000 <- get.sig.thr( perms_Log10wkACR_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr_Log15wkACR_1000 <- get.sig.thr( perms_Log15wkACR_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

#plot with significance threshold

tiff("./Plots/Col4a5xDO 6wk ACR QTL perms_1000.tiff")
plot(ACR_6wks_qtl, sig.thr = thr_6wkACR_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 6wk ACR QTL perms_1000")
dev.off()
tiff("./Plots/Col4a5xDO 10k ACR QTL perms_1000.tiff")
plot(ACR_10wks_qtl, sig.thr = thr_10wkACR_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 10wk ACR QTL perms_1000")
dev.off()
tiff("./Plots/Col4a5xDO 15wk ACR QTL perms_1000.tiff")
plot(ACR_15wks_qtl, sig.thr = thr_15wkACR_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 15wk ACR QTL perms_1000")
dev.off()

tiff("./Plots/Col4a5xDO 6wk Log[ACR] QTL perms_1000.tiff")
plot(Log_ACR_6wks_qtl, sig.thr = thr_Log6wkACR_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 6wk Log[ACR] QTL perms_1000")
dev.off()
tiff("./Plots/Col4a5xDO 10wk Log[ACR] QTL perms_1000.tiff")
plot(Log_ACR_10wks_qtl, sig.thr = thr_Log10wkACR_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 10wk Log[ACR] QTL perms_1000")
dev.off()
tiff("./Plots/Col4a5xDO 15wk Log[ACR] QTL perms_1000.tiff")
plot(Log_ACR_15wks_qtl, sig.thr = thr_Log15wkACR_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO 15wk Log[ACR] QTL perms_1000")
dev.off()

# Plotting coefficents by Chr ----------------------------------------
#Chr close up coefficient plot by CC founders

# ACR 6wk
ACR_6wks_qtl$coef$A[abs(ACR_6wks_qtl$coef$A) > 1000 ] = 0 
coefplot( GFR_qtl, chr = 2, main = "Col4a5xDO GFR QTL Chr2")
# ACR 10wk
temp_10wkACR <- ACR_10wks_qtl
ACR_10wks_qtl$coef$A[abs(ACR_10wks_qtl$coef$A) > 10000 ] = 0 
coefplot( ACR_10wks_qtl, chr = 13, main = "Col4a5xDO 10wk ACR QTL Chr13")
# ACR 15wk
temp_15wkACR <- ACR_15wks_qtl
ACR_15wks_qtl$coef$A[abs(ACR_15wks_qtl$coef$A) > 4000 ] = 0 
coefplot( ACR_15wks_qtl, chr = 2, main = "Col4a5xDO 15wk ACR QTL Chr2")

# Log[ACR] 6wk
temp_log_6wkACR <- Log_ACR_6wks_qtl
Log_ACR_6wks_qtl$coef$A[abs(Log_ACR_6wks_qtl$coef$A) > 100 ] = 0 
coefplot( Log_ACR_6wks_qtl, chr = 11, main = "Col4a5xDO 6wk Log[ACR] QTL Chr11")
# Log[ACR] 10wk
temp_log_10wkACR <- Log_ACR_10wks_qtl
Log_ACR_10wks_qtl$coef$A[abs(Log_ACR_10wks_qtl$coef$A) > 2 ] = 0 
coefplot( Log_ACR_10wks_qtl, chr = 16, main = "Col4a5xDO 10wk Log[ACR] QTL Chr16")
# Log[ACR] 15wk
temp_log_15wkACR <- Log_ACR_15wks_qtl
Log_ACR_15wks_qtl$coef$A[abs(Log_ACR_15wks_qtl$coef$A) > 2 ] = 0 
coefplot( Log_ACR_15wks_qtl, chr = 15, main = "Col4a5xDO 15wk Log[ACR] QTL Chr15")

# Interval peaks  -----------------

interval = bayesint(GFR_log_qtl, chr = )
knitr::kable(interval)

pxg.plot(pheno = pheno, pheno.col = "C2_log", probs = probs, snp.id = interval[2,1], snps = GM_snps)

# GWAS--------------------------------------------------------------------------------

Sex_covar <- as.matrix(Sex_covar)

ACR_6wk_GWAS = scanone.assoc(pheno = pheno, pheno.col = "ACR6WK", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
                     sdp.file = "./DO_Sanger_SDPs.txt.bgz", ncl = 4 )
ACR_10wk_GWAS = scanone.assoc(pheno = pheno, pheno.col = "ACR10WK", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
                             sdp.file = "./DO_Sanger_SDPs.txt.bgz", ncl = 4 )
ACR_15wk_GWAS = scanone.assoc(pheno = pheno, pheno.col = "ACR15WK", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
                             sdp.file = "./DO_Sanger_SDPs.txt.bgz", ncl = 4 )

LogACR_6wk_GWAS = scanone.assoc(pheno = pheno, pheno.col = "LogACR6WK", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
                             sdp.file = "./DO_Sanger_SDPs.txt.bgz", ncl = 4 )
LogACR_10wk_GWAS = scanone.assoc(pheno = pheno, pheno.col = "LogACR10WK", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
                             sdp.file = "./DO_Sanger_SDPs.txt.bgz", ncl = 4 )
LogACR_15wk_GWAS = scanone.assoc(pheno = pheno, pheno.col = "LogACR15WK", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
                             sdp.file = "./DO_Sanger_SDPs.txt.bgz", ncl = 4 )


plot(gwas, main = "GWAS Log GFR")

#Gwas_perms = assoc.map.perms(pheno = pheno, pheno.col = "C2_log", probs = probs,  addcovar = Sex_covar, snps = GM_snps, 
#                             snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz", nperm = 100)

Gwas_perms_1000 = assoc.map.perms(pheno = pheno, pheno.col = "C2_log", probs = probs,  addcovar = Sex_covar, snps = GM_snps, 
                                  snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz", nperm = 1000)

#Gwas_perms_100 <- Gwas_perms

#Gwas_100_thr = get.sig.thr(-log10(Gwas_perms_100), Xchr = FALSE)
Gwas_1000_thr = get.sig.thr(-log10(Gwas_perms_1000), Xchr = FALSE)

#For 100 perms
#plot(gwas, sig.thr = Gwas_100_thr, main = "GWAS Log GFR perm100")
#plot(gwas, chr = 4, sig.thr = Gwas_100_thr,main = "GWAS Log GFR Chr4")
#plot(gwas, chr = 7, sig.thr = Gwas_100_thr,main = "GWAS Log GFR Chr7")
#plot(gwas, chr = 19, sig.thr = Gwas_100_thr, bin.size = 1, main = "GWAS Log GFR Chr19")

#For 1000 perms
plot(gwas, sig.thr = Gwas_1000_thr, main = "GWAS Log GFR perm1000")



# Candidate Genes---------------------------------------------------------------------------
# Chr7====
chr = 7
assoc = assoc.map(pheno = pheno, pheno.col ="C2_log", probs = probs, K = K_GM[[chr]],
                  addcovar = Sex_covar, snps = GM_snps, chr = chr, start = 121.2 ,
                  end = 122.5 , output = "p-value")

chr7_gene <- assoc

Chr7_gene_plot = assoc.plot(chr7_gene, thr = 5, show.sdps = FALSE)

# Chr 4 ======
chr = 4
assoc = assoc.map(pheno = pheno, pheno.col ="C2_log", probs = probs, K = K_GM[[chr]],
                  addcovar = Sex_covar, snps = GM_snps, chr = chr, start = 7 ,
                  end = 13 , output = "p-value")

chr4_gene <- assoc

Chr4_gene_plot = assoc.plot(chr4_gene, thr = 4, show.sdps = TRUE)

# Chr 19 =======
chr = 19
assoc = assoc.map(pheno = pheno, pheno.col ="C2_log", probs = probs, K = K_GM[[chr]],
                  addcovar = Sex_covar, snps = GM_snps, chr = chr, start =26 ,
                  end = 29 , output = "p-value")

chr19_gene <- assoc

Chr4_gene_plot = assoc.plot(chr19_gene, thr = 5, show.sdps = TRUE)
