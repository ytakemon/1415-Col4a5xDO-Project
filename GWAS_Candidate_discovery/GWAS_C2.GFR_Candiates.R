#GWAS and candidate genes under Chr 2,4,7, 9 and 19
library(DOQTL)
library(abind)
library(knitr)

setwd("/hpcdata/ytakemon/Col4a5xDO")

#load sample
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#read and clean up phenotype data
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of GFR
pheno[pheno < 0 ] = NA
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2) 
pheno$X2xC1_log <- log(pheno$X2xC1)
pheno$C1_log <- log(pheno$C1)
pheno$PL_log <- log(pheno$PL)
options(na.action = 'na.pass')
#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"

#create a subset separating female and male samples
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

#kinship mapping
K_GS <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = TRUE)
K_GS.F <- kinship.probs(genoprob.F, snps = GM_snps, bychr = TRUE)
K_GS.M <- kinship.probs(genoprob.M, snps = GM_snps, bychr = TRUE)

save(K_GS, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
save(K_GS.F, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.F.Rdata")
save(K_GS.M, file ="./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.M.Rdata")


#GWAS
#males <- which(sex.covar[,1] == 1)
#females <- which(sex.covar[,1] == 0)

#males.perm <- sample( rownames( pheno)[males])
#females.perm <- sample( rownames(pheno)[females])

#rownames(pheno)[males] <- males.perm 
#rownames(pheno)[females] <- females.perm 

#rownames(sex.covar)[males] <- males.perm
#rownames(sex.covar)[females] <- females.perm

sex.covar <- as.matrix(sex.covar)
GWAS.log.C2 <- scanone.assoc(
	pheno = pheno, 
	pheno.col = "C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS, 
	addcovar = sex.covar, 
	markers = GM_snps, 
	sdp.file = "./GBRS_reconstruction/reconstruct/resources/DO_Sanger_SDPs.txt.bgz", 
	ncl = 6)
save(GWAS.log.C2, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas/gwas.log.C2.Rdata")

#creating a permutation then run this 1000 times
#min.p <- min(sapply(GWAS.log.C2[1:19], function(z){min(z$p.value)}))
#write(min.p, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm/perm1.txt")
#jump to unix
#grep . *.txt > all.perm.txt

#back to R
GFR.C2.perm <- read.delim("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/gwas.perm/gfr.c2.perm/all.perm.txt", sep ="\t", header = F )
gwas.thr.1000 <- get.sig.thr(-log10(GFR.C2.perm), alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.log.C2.GFR.png", width = 2000, height = 1000, res = 100)
plot(GWAS.log.C2, ylim = c(0 ,max(gwas.thr.1000)), sig.thr = gwas.thr.1000, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO GWAS of log C2 GFR")
dev.off()

#
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/GWAS.chr.log.C2.GFR.png", width = 1500, height = 1000, res = 100)
layout(matrix(1:8,2,4))
plot(GWAS.log.C2 , chr = 2, main = "Chr 2 log.C2.GFR")
plot(GWAS.log.C2 , chr = 4, main = "Chr 4 log.C2.GFR")
plot(GWAS.log.C2 , chr = 7, main = "Chr 7 log.C2.GFR")
plot(GWAS.log.C2 , chr = 9, main = "Chr 9 log.C2.GFR")
plot(GWAS.log.C2 , chr = 10, main = "Chr 10 log.C2.GFR")
plot(GWAS.log.C2 , chr = 15, main = "Chr 10 log.C2.GFR")
plot(GWAS.log.C2 , chr = 19, main = "Chr 19 log.C2.GFR")
dev.off()


######################################
#Create allele effect map for Chr 2,4,7,9,10,15,19
#coefficient threshold
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata" ) # qtl.GFR.log.C2.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata") # K_GS
qtl <- qtl.GFR.log.C2.192

qtl$coef$A[abs(qtl$coef$A) > 5 ] = 0 
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 2, main = "Chr 2 log.C2.GFR")
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 10 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr4.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 4, main = "Chr 4 log.C2.GFR")
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 10 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr7.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 7, main = "Chr 7 log.C2.GFR")
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 1 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr9.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 9, main = "Chr 9 log.C2.GFR")
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 1000 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr10.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 10, main = "Chr 10 log.C2.GFR")
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 1000 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr15.check.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 15, main = "Chr 15 log.C2.GFR")
dev.off()
qtl$coef$A[abs(qtl$coef$A) > 1000 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.log.C2.GFR.chr19.png", width = 1500, height = 1000, res = 100)
coefplot(qtl , chr = 19, main = "Chr 19 log.C2.GFR")
dev.off()

###################################
#Candidate gene hunting....
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata" ) # qtl.GFR.log.C2.192
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata") # K_GS
qtl <- qtl.GFR.log.C2.192
#create interval
interval = bayesint(qtl, chr = 2)
knitr::kable(interval)

#Hunting.......
#chr2
#|            |marker      |chr |        pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|----------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|UNCHS003889 |UNCHS003889 |2   |   5.608243|  1.261172| 10.840040| 14.22753| 3.089468| 0.0472794|   1.3253279|
#|UNC3158171  |UNC3158171  |2   |  60.316984| 29.708598| 17.250941| 23.48034| 5.098690| 0.0014053|   2.8522156|
#|UNC4291363  |UNC4291363  |2   | 156.601026| 69.896277|  9.110123| 11.84467| 2.572038| 0.1057840|   0.9755801|
interval = bayesint(qtl, chr = 2)
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = chr, 
	start = 20,
	end = 40, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr2.narrow.20.40.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 10, show.sdps = TRUE)
dev.off()

interval = bayesint(qtl, chr = 2)
chr <- 2
chr2.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = chr, 
	start = 55,
	end = 65, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr2.narrow.55.65.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr2.genes, thr = 10, show.sdps = TRUE)
dev.off()

#chr4
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|JAX00115386 |JAX00115386 |4   | 3.209815| 0.019268| 19.36596| 26.69091| 5.795858| 0.0003788|    3.421622|
#|JAX00115682 |JAX00115682 |4   | 7.187726| 0.892403| 20.81974| 28.94695| 6.285751| 0.0001479|    3.829900|
#|UNC6691401  |UNC6691401  |4   | 8.203528| 1.351902| 15.16022| 20.38630| 4.426830| 0.0047930|    2.319391|

interval = bayesint(qtl, chr = 4)
chr <- 4
chr4.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = chr, 
	start = 0,
	end = 10, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr4.narrow.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr4.genes, thr = 4, show.sdps = TRUE)
dev.off()

#chr7
#|            |marker      |chr |       pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|UNC12559491 |UNC12559491 |7   |  29.60412|  9.178563|  8.696113| 11.28113| 2.449665| 0.1268179|   0.8968196|
#|JAX00653739 |JAX00653739 |7   | 121.88554| 55.308275| 18.655894| 25.60374| 5.559782| 0.0005927|   3.2271615|
#|UNC13969401 |UNC13969401 |7   | 138.36300| 72.967777|  8.786200| 11.40353| 2.476246| 0.1219600|   0.9137825|

interval = bayesint(qtl, chr = 7)
chr <- 7
chr7.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 7, 
	start = 120,
	end = 126, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr7.candiates.narrow.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr7.genes, thr = 1.85, show.sdps = TRUE)
dev.off()

#chr9
#|            |marker      |chr |      pos|        cM|  perc.var|       lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|---------:|---------:|---------:|--------:|---------:|-----------:|
#|UNC16083432 |UNC16083432 |9   | 30.33947|  9.424627|  8.751890| 11.356901| 2.466120| 0.1237911|   0.9073106|
#|UNCHS025062 |UNCHS025062 |9   | 32.61934| 12.172372| 16.873873| 22.916580| 4.976272| 0.0017626|   2.7538541|
#|UNC16522104 |UNC16522104 |9   | 65.35888| 32.188614|  7.689356|  9.921331| 2.154390| 0.1930756|   0.7142726|

interval = bayesint(qtl, chr = 9)
chr <- 9
chr9.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 9, 
	start = 30,
	end = 40, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr9.narrow.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr9.genes, thr = 1.85, show.sdps = TRUE)
dev.off()

#chr10
#|            |marker      |chr |        pos|        cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|----------:|---------:|---------:|--------:|--------:|---------:|-----------:|
#|ICR4223     |ICR4223     |10  |   7.103149|  0.051242|  8.023105| 10.37046| 2.251918| 0.1685334|   0.7733141|
#|UNC18427765 |UNC18427765 |10  |  91.435548| 39.127881| 16.938049| 23.01235| 4.997068| 0.0016962|   2.7705261|
#|JAX00302382 |JAX00302382 |10  | 128.348227| 65.644214|  7.858797| 10.14915| 2.203860| 0.1802777|   0.7440579|

interval = bayesint(qtl, chr = 10)
chr <- 10
chr10.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 10, 
	start = 88,
	end = 98, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr10.narrow.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr10.genes, thr = 3, show.sdps = TRUE)
dev.off()

interval = bayesint(qtl, chr = 15)
knitr::kable(interval)

#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNC25475383 |UNC25475383 |15  |  43.98596| 13.89354|  9.66417| 12.60287| 2.736677| 0.0823963|    1.084092|
#|UNC25582381 |UNC25582381 |15  |  51.66936| 16.94894| 18.05776| 24.69529| 5.362513| 0.0008592|    3.065923|
#|UNC26247792 |UNC26247792 |15  | 103.65408| 54.49084| 13.46949| 17.93947| 3.895506| 0.0122461|    1.912002|
#pt1
chr <- 15
chr15.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 15, 
	start = 50,
	end = 55, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/gfr.chr15.pt1.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr15.genes, thr = 10, show.sdps = TRUE)
dev.off()
#pt2
chr <- 15
chr15.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 15, 
	start = 85,
	end = 90, 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/gfr.chr15.pt2.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr15.genes, thr = 10, show.sdps = TRUE)
dev.off()

#chr15 
#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNC25475383 |UNC25475383 |15  |  43.98596| 13.89354|  9.66417| 12.60287| 2.736677| 0.0823963|    1.084092|
#|UNC25582381 |UNC25582381 |15  |  51.66936| 16.94894| 18.05776| 24.69529| 5.362513| 0.0008592|    3.065923|
#|UNC26247792 |UNC26247792 |15  | 103.65408| 54.49084| 13.46949| 17.93947| 3.895506| 0.0122461|    1.912002|

interval = bayesint(qtl, chr = 19)
chr <- 19
chr19.genes <- assoc.map(
	pheno = pheno, 
	pheno.col ="C2_log", 
	probs = best.genoprobs.192, 
	K = K_GS[[chr]],
	addcovar = sex.covar, 
	snps = GM_snps, 
	chr = 19, 
	start = interval[1,3],
	end = interval[3,3], 
	output = "p-value")
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Chr19.candiates.png", width = 1500, height = 1000, res = 100)
assoc.plot(chr19.genes, thr = 5, show.sdps = TRUE)
dev.off()


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
##############################################################Check Genotype Distribution##################################################################

#Chr15
interval = bayesint(qtl, chr = 15)
knitr::kable(interval)
#|            |marker      |chr |       pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNC25475383 |UNC25475383 |15  |  43.98596| 13.89354|  9.66417| 12.60287| 2.736677| 0.0823963|    1.084092|
#|UNC25582381 |UNC25582381 |15  |  51.66936| 16.94894| 18.05776| 24.69529| 5.362513| 0.0008592|    3.065923|
#|UNC26247792 |UNC26247792 |15  | 103.65408| 54.49084| 13.46949| 17.93947| 3.895506| 0.0122461|    1.912002|

GM_15 <- GM_snps[GM_snps$chr == 15,]
GM_15.89 <- GM_15[GM_15$pos < 88.2,]
GM_15.87 <- GM_15.89[GM_15.89$pos > 87.7,]


#look at pos = 87.73390 UNC26053317
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/gfr.chr15.geno.dist.87.7.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, snp.id = "UNC26053317", snps = GM_snps)
dev.off()
#look at post = 88.18756 JAX00408903
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/gfr.chr15.geno.dist.88.1.png", width = 1500, height = 1000, res = 100)
pxg.plot(pheno = pheno, pheno.col = "C2_log", probs = best.genoprobs.192, snp.id = "JAX00408903", snps = GM_snps)
dev.off()






















