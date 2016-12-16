#source and laod DOQTL package: this only needs to happen once
#source("http://bioconductor.org/biocLite.R")
#biocLite("DOQTL")
#library(devtools)
#library(knitr)
#install_github("dmgatti/DOQTL", force = TRUE)


library(DOQTL)
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/")

# Load data -------------------------------------------

#load geno probs data
load("~/Desktop/Col4a5xDO/R_col4a5xdo/new_col4a5_probs.Rdata")

#create kinship probablility from JAX ftp site containg GigaMuga array snps. 
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#remove I for 9th genotype, as this will mess up the 8 founder IDs
colnames(probs) = sub("I", "", colnames(probs)) 

# Clean up data -----------------------------------------------

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
probs = probs[-c(31,75,80,81,84,85,86),,] #end result 175 samples left

#load in phenotype data containing sex, generation, GFR, and ACR
pheno <- read.table("../Phenotype/1415_master_pheno.txt", header = TRUE, sep= "\t")
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names

#GFR C2 has negative values so change to NA
pheno[pheno < 0 ] = NA
pheno$C2[123] = NA #over sigma 650000 cut off
pheno$C2_log = log(pheno$C2)


#extract samples that exist in both pheno and probs should be 175 samples
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

# Not sure what this is at the momement ---------------------
###################
#duplicates
tiff("./Plots/Duplicated Samples.tiff")
layout(matrix(1:6,1,2))
image(probs[24,,], main = rownames(probs)[24])
image(probs[31,,], main = rownames(probs)[31])
dev.off()

#B6 like samples
tiff("./Plots/B6-like samples.tiff")
layout(matrix(1:6,2,3))
image(probs[75,,], main = rownames(probs)[75])
image(probs[80,,], main = rownames(probs)[80])
image(probs[81,,], main = rownames(probs)[81])
image(probs[84,,], main = rownames(probs)[84])
image(probs[85,,], main = rownames(probs)[85])
image(probs[86,,], main = rownames(probs)[86])
dev.off()

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

#Scanone for GFR qtl
GFR_qtl = scanone( pheno = pheno, pheno.col = "C2_log", probs = probs, K = K_GM,
                   addcovar = Sex_covar, snps = GM_snps)
GFR_log_qtl <- GFR_qtl
#GFR_temp = GFR_qtl

#plot qtl
plot(GFR_qtl, main ="Col4a5xDO GFR_log QTL")

# Create permuations for significante for plotting -------------------------------

#cerate permuations for significance testing
perms_log_100 = scanone.perm(pheno = pheno, pheno.col = "C2_log", probs = probs,  
                          addcovar = Sex_covar, snps = GM_snps, nperm = 100)

perms_log_1000 = scanone.perm(pheno = pheno, pheno.col = "C2_log", probs = probs,  
                             addcovar = Sex_covar, snps = GM_snps, nperm = 1000)


#perms_1000 = scanone.perm(pheno = pheno, pheno.col = "C2_", probs = probs,  
#                     addcovar = Sex_covar, snps = GM_snps, nperm = 1000)

#perms_1000 = scanone.perm(pheno = pheno, pheno.col = "GFR_C2", probs = model.probs, addcovar = Sex_covar,
#                     snps = MM_snps, nperm = 1000)

#create threshold for alphas
thr_log_100 = get.sig.thr( perms_log_100[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

thr_log_1000 = get.sig.thr( perms_log_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)


#thr = get.sig.thr( perms_100[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE) #get.sig.thr not working
#thr_100 <- thr
#thr_1000 = get.sig.thr( perms_1000[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE) #get.sig.thr not working

#raw GFRC2 data, not log normalized
#plot(GFR_qtl, sig.thr = thr_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GFR QTL perms_1000")
#plot(GFR_qtl, sig.thr = thr_100, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GFR QTL perms_100")


#log normalized data
plot(GFR_qtl, sig.thr = thr_log_100, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GFR_log QTL perms_100")
plot(GFR_qtl, sig.thr = thr_log_1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GFR_log QTL perms_1000")

# Plotting coefficents by Chr ----------------------------------------

#coefficient threshold
GFR_qtl$coef$A[abs(GFR_qtl$coef$A) > 1000 ] = 0 

#Chr close up coefficient plot by CC founders
# plot for raw GFR C2 data
coefplot( GFR_qtl, chr = 2, main = "Col4a5xDO GFR QTL Chr2")
coefplot( GFR_qtl, chr = 3, main = "Col4a5xDO GFR QTL Chr3")
coefplot( GFR_qtl, chr = 6, main = "Col4a5xDO GFR QTL Chr6")
coefplot( GFR_qtl, chr = 12, main = "Colx4a5xDO GFR QTL Chr12")
coefplot( GFR_qtl, chr = 15, main = "Colx4a5xDO GFR QTL Chr15")
# plot for log transformed GFR C2 data
coefplot( GFR_qtl, chr = 2, main = "Col4a5xDO Log GFR QTL Chr2")
coefplot( GFR_qtl, chr = 4, main = "Col4a5xDO Log GFR QTL Chr4")
coefplot( GFR_qtl, chr = 7, main = "Col4a5xDO Log GFR QTL Chr7")
coefplot( GFR_qtl, chr = 15, main = "Col4a5xDO Log GFR QTL Chr15")
coefplot( GFR_qtl, chr = 19, main = "Col4a5xDO Log GFR QTL Chr19")

# Interval peaks  ---------------------------------------------------------

interval = bayesint(GFR_log_qtl, chr = 19)
knitr::kable(interval)

pxg.plot(pheno = pheno, pheno.col = "C2_log", probs = probs, snp.id = interval[2,1], snps = GM_snps)

# GWAS--------------------------------------------------------------------------------

Sex_covar <- as.matrix(Sex_covar)

gwas = scanone.assoc(pheno = pheno, pheno.col = "C2_log", probs = probs, K = K_GM, addcovar = Sex_covar, markers = GM_snps, 
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
