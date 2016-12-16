## Col4a5xDO eQTL analysis
#	Author: Yuka Takemon
#	Date Created: 11/30/16
#	Date Last Modified: 12/06/16

## Ron identified genes to follow up with for eQTL analysis, so that we can narrow down to a few genes to send to Dan our new
#	collaborator who has a human cohort of Alport patients from Cyprus.
#	We will do doign an eQTL analysis of the following genes:
## From the GFR qtl:
#	Kcnv2 	(chr19)	ENSMUSG00000047298
#	Pum3 	(chr19)	ENSMUSG00000041360
#	Rfx3 	(chr19)	ENSMUSG00000040929
## From the Alb/ACR 6wk :
#	Dgke	(Chr11)	ENSMUSG00000000276
## From the Alb/ACR 10wk :
#	Pik3r1	(Chr13)	ENSMUSG00000041417
## From the Alb/ACR 6 and 10:
#	Fmn1 		(chr2)	ENSMUSG00000044042
## From the ACR 15wk:
#	Tbc1d22a	(chr15)	ENSMUSG00000051864
#	Fam19a5		(chr15)	ENSMUSG00000054863
#	Brd1		(chr15) ENSMUSG00000022387
#	Zbed4		(chr15)	ENSMUSG00000034333
#	Alg12		(chr15)	ENSMUSG00000035845
#	Creld2		(chr15)	ENSMUSG00000023272
library(DOQTL)
library(knitr)
setwd("/hpcdata/ytakemon/Col4a5xDO")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/K_GS.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples

#create sex covariate
sex.covar <- model.matrix(~0+Sex, data = pheno)
colnames(sex.covar)[2] <- "sex"
sex.covar <- sex.covar[,"sex"]
sex.covar <- as.data.frame(sex.covar)
colnames(sex.covar)[1] <- "sex"
save(sex.covar, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/sex.covar.Rdata")

#RankZ
RNA_seqZ <- t(RNA_seq)
RNA_seqZ <- apply(RNA_seqZ, 2, rankZ)
RNA_seqZ <- t(RNA_seqZ)

##GFR Chr19 candidates
#	QTL Kcnv2
qtl.Kcnv2 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000047298", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Pum3
qtl.Pum3	 	<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000041360", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Rfx3
qtl.Rfx3 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000040929", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)

##Alb 6wk Chr 11 genenes
#	QTL Dgke
qtl.Dgke 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000000276", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)

##Alb110wk Chr 15 genes
#	QTL Pik3r1
qtl.Pik3r1		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000041417", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)

##ACR 15wk Chr15 genes
#	QTL Tbc1d22a
qtl.Tbc1d22a 	<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000051864", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Fam19a5
qtl.Fam19a5 	<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000054863", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Brd1
qtl.Brd1		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000022387", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Zbed4
qtl.Zbed4 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000034333", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Alg12
qtl.Alg12 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000035845", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)
#	QTL Creld2
qtl.Creld2 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000023272", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)

##ACR 10wk Chr2 gene
#	QTL Fmn1 	
qtl.Fmn1 		<- scanone( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000044042", probs = best.genoprobs.192, K = K_GS, addcovar = sex.covar, snps = GM_snps)


save(qtl.Kcnv2, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Kcnv2.Rdata")
save(qtl.Pum3, 		file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Pum3.Rdata")
save(qtl.Rfx3, 		file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Rfx3.Rdata")
save(qtl.Tbc1d22a, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Tbc1d22a.Rdata")
save(qtl.Fam19a5, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fam19a5.Rdata")
save(qtl.Brd1,		file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Brd1.Rdata")
save(qtl.Zbed4, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Zbed4.Rdata")
save(qtl.Alg12, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Alg12.Rdata")
save(qtl.Creld2, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Creld2.Rdata")
save(qtl.Fmn1, 		file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1.Rdata")
save(qtl.Dgke, 		file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Dgke.Rdata")
save(qtl.Pik3r1, 	file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Pik3r1.Rdata")

##PERMUTATION
# This process needs to only happen once with one gene, and the can be applied to all genes.
#perms.1000.qtl.RNA.rankZ.tpm  <- scanone.perm( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000047298", probs = best.genoprobs.192, addcovar = sex.covar, snps = GM_snps, nperm = 1000)
#save(perms.1000.qtl.RNA.rankZ.tpm, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.RNA.rankZ.tpm.Rdata" )
#load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.RNA.rankZ.tpm.Rdata")
#threshold <- get.sig.thr( perms.1000.qtl.RNA.rankZ.tpm[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)
#save(threshold, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Kcnv2.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Pum3.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Rfx3.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Tbc1d22a.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fam19a5.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Brd1.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Zbed4.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Alg12.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Creld2.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Dgke.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Pik3r1.Rdata")

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Kcnv2.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Kcnv2, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Kcnv2 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Pum3.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Pum3, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Pum3 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Rfx3.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Rfx3, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Rfx3 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Tbc1d22a.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Tbc1d22a, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Tbc1d22a QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fam19a5.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Fam19a5, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Fam19a5 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Brd1.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Brd1, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Brd1 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Zbed4.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Zbed4, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Zbed4 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Alg12.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Alg12, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Alg12 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Creld2.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Creld2, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Creld2 QTL perm 1000")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1.qtl.perm1000.png", width = 1500, height = 1000, res = 100)
plot(qtl.Fmn1, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Fmn1 QTL perm 1000")
dev.off()
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Dgke.qtl.perm1000.pdf", width = 10.0, height = 7.5)
plot(qtl.Dgke, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Dgke QTL perm 1000")
dev.off()
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Pik3r1.qtl.perm1000.pdf", width = 10.0, height = 7.5)
plot(qtl.Pik3r1, sig.thr = threshold, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO Pik3r1 QTL perm 1000")
dev.off()

######################


## Post meeting with Ron we decided to pursue Rfx3 (Chr 19 from GFR QTL) and Fmn1 (Chr 2 from ACR/Alb at 6 and 10 wk QTL)
#################################################################
#	Coefficient Allele effect plot
#	load qtl data 
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Rfx3.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/qtl.RNA.rankZ.tpm.Fmn1.Rdata")

#	Allele effect plots
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Rfx3.coef.chr19.png", width = 1500, height = 1000, res = 100)
coefplot(qtl.Rfx3, chr = 19, main = "Allele effect of Rfx3 QTL at Chr 19")
dev.off()

qtl.Fmn1$coef$A[abs(qtl.Fmn1$coef$A) > 0.05 ] = 0
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1.coef.chr2.png", width = 1500, height = 1000, res = 100)
coefplot(qtl.Fmn1, chr = 2, main = "Allele effect of Fmn1 QTL at Chr 2")
dev.off()

#	Allele distribution of region of interst
#	load qtl data from previous analysis that had the narrowest bayesian interval to help narrow down region of interest
#	Rfx3 		ENSMUSG00000040929
#	Fmn1 		ENSMUSG00000044042
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")

#Rfx3
GFR_interval 	<- bayesint(qtl.GFR.log.C2.192, chr = 19)
Rfx3_interval 	<- bayesint(qtl.Rfx3, chr = 19)
knitr::kable(GFR_interval)
knitr::kable(Rfx3_interval)
#> knitr::kable(GFR_interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNC30136851 |UNC30136851 |19  | 26.98658| 14.99423| 20.10076| 27.82608| 6.042356| 0.0002364|    3.626293|
#|JAX00474835 |JAX00474835 |19  | 27.61228| 15.38569| 24.29844| 34.51806| 7.495501| 0.0000138|    4.860904|
#|UNCHS047487 |UNCHS047487 |19  | 28.04451| 15.68458| 16.49212| 22.34842| 4.852898| 0.0022117|    2.655270|
#> knitr::kable(Rfx3_interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNC30135096 |UNC30135096 |19  | 26.87585| 14.79494| 12.52813| 25.69975| 5.580631| 0.0005698|    3.244270|
#|UNCHS047465 |UNCHS047465 |19  | 27.60969| 15.37632| 16.80785| 35.33131| 7.672096| 0.0000097|    5.013693|
#|UNCHS047480 |UNCHS047480 |19  | 27.91010| 15.53331| 14.89991| 30.97768| 6.726718| 0.0000628|    4.202276|
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Rfx3.allele.dist1.png", width = 1500, height = 1000, res = 100)
pxg.plot( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000040929", probs = best.genoprobs.192, snp.id = "JAX00474835", snps = GM_snps)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Rfx3.allele.dist2.png", width = 1500, height = 1000, res = 100)
pxg.plot( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000040929", probs = best.genoprobs.192, snp.id = "UNCHS047465", snps = GM_snps)
dev.off()

#Fmn1
Alb_interval	<- bayesint(qtl.log.Alb10WK.192, chr = 2)
Fmn1_interval	<- bayesint(qtl.Rfx3, chr = 2)
knitr::kable(Alb_interval)
knitr::kable(Fmn1_interval)
#> knitr::kable(Alb_interval)
#|            |marker      |chr |      pos|       cM| perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|--------:|--------:|--------:|--------:|--------:|---------:|-----------:|
#|UNCHS005963 |UNCHS005963 |2   | 101.7360| 45.94599| 11.89948| 20.39737| 4.429233| 0.0047724|    2.321265|
#|JAX00098817 |JAX00098817 |2   | 113.4217| 50.67787| 16.48704| 29.00710| 6.298813| 0.0001443|    3.840867|
#|UNC3776734  |UNC3776734  |2   | 113.6730| 50.81143| 11.95717| 20.50284| 4.452134| 0.0045801|    2.339128|
#> knitr::kable(Fmn1_interval)
#|            |marker      |chr |       pos|       cM|  perc.var|      lrs|      lod|         p| neg.log10.p|
#|:-----------|:-----------|:---|---------:|--------:|---------:|--------:|--------:|---------:|-----------:|
#|UNCHS004963 |UNCHS004963 |2   |  52.65429| 26.22784|  7.826348| 15.64721| 3.397747| 0.0285424|    1.544509|
#|UNCHS004986 |UNCHS004986 |2   |  54.74089| 26.79108| 11.684956| 23.85787| 5.180670| 0.0012068|    2.918375|
#|UNC4519006  |UNC4519006  |2   | 172.11418| 84.14639|  8.591619| 17.24794| 3.745342| 0.0158655|    1.799546|
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1.allele.dist1.png", width = 1500, height = 1000, res = 100)
pxg.plot( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000044042", probs = best.genoprobs.192, snp.id = "JAX00098817", snps = GM_snps)
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Fmn1.allele.dist2.png", width = 1500, height = 1000, res = 100)
pxg.plot( pheno = RNA_seqZ, pheno.col = "ENSMUSG00000044042", probs = best.genoprobs.192, snp.id = "UNC3776734", snps = GM_snps)
dev.off()

















