##Xce and Xist allele for Col4a5xDO
#	Yuka Takemon
#	Created: 11/09/16
#	Modified last: 12/07/16

##Intro:
#	We know that X chromosome silencing is driven by a region called the Xce (X controlling element), which is a cis-regulatory gene.
#	Xist, 5' to Xce, is expressed by the inactive X-chromosome and facilitates the spread of inactivation of that chromosome.
#	In Mus musculus there are 4 known Xce alleles: Xce^a, Xce^b, Xce^c, and Xce^e, which effect the probability that the X-chromosome will be
#	expressed. The liklihood is as follows: Xce^a < Xce^e < Xce^b < Xce^c
#	We would like to perform ANOVA test to see if the allele of Xce allele has an effect on the phenotype on females as they are Hets for the Col4a5
#	mutation located on the B6 X-chromosome.

##Xce allels of 8 founders:
#	Xce^a : AJ, 129
#	Xce^b : WSB, B6, NOD, NZO
#	Xce^c : CAST
#	Xce^e :	PWK

##Expectation of phenotype severity:
#	Xce^a x Xce^b: will have more severe phenotype as the mutation is on the Xce^b chromosome.
#	Xce^b x Xce^b: will have equal distribution of phenotype as the probability is equal.
#	Xce^c x Xce^b: will have less severe phenotype as the dominant allele does not carry the mutaiton.
#	Xce^e x Xce^b: will have more severe phenotype as the Xce^e is less dominanat that the allele carrying the mutation.

##Xce is not yet on Ensemble but there are a few papers that have located a specific locus that suggests to the the Xce alelle. (Calaway et al., 2013)
#	Using the two markers: JAX00716733 and JAX00182728 we can distinguish between the four alleles.
#	Xce allele JAX00716733 JAX00182728
#	Xce^a 		C 			C
#	Xce^b 		C 			T
#	Xce^c 		G 			G
#	Xce^e 		G 			T

##Xist vs Xce
#	We can check the effects of Xce by identifying the gene expression of Xist gene and which gene it was derived from as the EMASE is a able to account
#	for this. Xist should be expressed form the X chromosoem that is being suppressed (the weaker Xce allele)

##Workflow:
#	First identify the snp region and assign specific alleles at the two markers using Genoprobs data (I think)
#	Then using the EMSE output identify founders that contirubted to the Xist expression
#	Compare Xce and Xist
#	Compare Xce to phenotypes (GFR & ACR)
#	Compare Xist to phenotypes (GFR & ACR)

##Xce interval
#	x:102826982-102990039

#setup:
library(DOQTL)
library(abind)
library(knitr)
library(ggplot2)
library(grid)

dir <- "/projects/marralab/ytakemon_prj/Col4a5/"
#load sample
load(paste0(dir, "Data/consolidated/best.genoprobs.192.Rdata"))
load(paste0(dir, "Data/consolidated/GM_snps.Rdata"))

##Identify Xce region and find GigaMuga array markers within this region:
#	In Calaway et al., 2013, the MDA snps show an interval of : x:102826982-102990039
#	We need to first find the snp markers in the GigaMuga array that are within this region.

#subset to only the X chromosome
GM_X <- GM_snps[GM_snps$chr == "X",]

#look up interval region
 x <- GM_X[GM_X$pos < 102.9990,]
 x <- x[x$pos > 102.8000,]

#	Bellow are markers within this interval
#XiD1, 102.8279 (T G)
#XiD2, 102.8283 (T G)
#XiE1 (C A)
#XiE2 (C A)
#UNC31159403, 102.8689, rs263068729 (T C) X:102868903
#JAX00716738, 102.8959, rs29082023 (C A)  X:102895879
#UNC31159628 102.9165 rs31398824 (G A) X:102916497
#UNC31159865 102.9735 rs264074348 (G A) X:102973501
#JAX00239728 102.9839 rs49189402 (C T) X:102983920

##pick 5 animals from the genoprob file and look up each of the 9 snp markers, theoretically they should all have the same genoprob call.
#	the following 5 animals are used for verification.
#"X1415.0082"
#"X1415.0134"
#"X1415.0540"
#"X1415.0578"
#"X1415.0864"

#best.genoprobs.192 was created by subseting dim to match Gillian's dataset. As a result some only 4 markers are in the best.genoprobs.192 dataset.
#The snp markers to look up are:
#UNC31159403
#JAX00716738
#UNC31159628
#UNC31159865

sample_names <- rownames(best.genoprobs.192)
marker_names <- c("UNC31159403","JAX00716738","UNC31159628","UNC31159865")

Xce_markers <- array(0, c(192, 2, 4), dimnames  = list( sample_names, c("max.founders", "max.probs"), marker_names)) #192 samples, 2col(allele, maxprob), 4 markers
for (sample in sample_names){

	for (snp in marker_names){

		df <- as.data.frame(best.genoprobs.192[sample,, snp])
		max.founder <- colnames(best.genoprobs.192)[apply(df,2,which.max)] #2 for column, #1 for row
		max.prob <- max(df)

		Xce_markers[sample, "max.founders", ] <- max.founder
		Xce_markers[sample,	"max.probs", ] <- max.prob

	}
}

rownames(Xce_markers)[Xce_markers[,"max.probs",] < 0.95] #samples less than 95% max.prob at Xce markers

samples_under <- c("X1415.0101", "X1415.0176", "X1415.0300")

for (sample in samples_under){
	print(sample)
	for (marker in marker_names){
		x <- best.genoprobs.192[sample, , marker]
		print(marker)
		print(x)
}
}

#	X1415.0101 is a B (0.2-0.3) & E (0.7-0.8) het, which is fine because both B6 and NZO have XceB allele.
#	X1415.0176 is also a B (0.2 - 0.4) & E (0.6 - 0.8) het.
#	X1415.0300 is a A (0.4) & C (0.6) het, AJ & 129 respectively, both of which have the XceA allele.
#	since they are all going to be called the same Xce allele I am moving on without correction

#	Assign Xce a, b, c, and e to each animal

for (sample in samples_under){
	x <- Xce_markers[sample,,]
	print(sample)
	print(x)
}

Xce_allele <- array(0, c(192,1), dimnames = list ( sample_names, "Xce_allele"))
for (sample in sample_names){
	founder <- unique(Xce_markers[ sample, "max.founders",])
	Xce_allele[sample, ] <- founder

	if ((founder == "A") || (founder == "C")){
		Xce_allele[sample, ] <- "Xce_a"
	}

	else if ((founder == "B") || (founder == "D") || (founder == "E") || (founder == "H")){
		Xce_allele[sample, ] <- "Xce_b"
	}

	else if ((founder == "F")){
		Xce_allele[sample, ] <- "Xce_c"
	}

	else if ((founder == "G")){
		Xce_allele[sample, ] <- "Xce_e"
	}
	else {
		print("error with sample: ", sample)
	}
}

Xce_allele <- as.data.frame(Xce_allele)
#This is the DO Xce allele, now accounting for B6 allele figure out which Xce is most likely to be expressed.
#Subset females out
pheno <- read.delim(paste0(dir,"Data/consolidated/Phenotype/1415_master_pheno.txt"), sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(Xce_allele),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$C2[76] = NA #over sigma 650000 cut off
pheno$C2[138] = NA #over sigma 650000 cut off

pheno$C2_log <- log(pheno$C2)

pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR10WK_log <- log(pheno$ACR10WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)

pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Alb10WK_log <- log(pheno$Alb10WK)
pheno$Alb15WK_log <- log(pheno$Alb15WK)

pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno$Creat10WK_log <- log(pheno$Creat10WK)
pheno$Creat15WK_log <- log(pheno$Creat15WK)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs

#subsetdat to only females
pheno.F <- subset(pheno, Sex == "F")
#Xce.data gives me the DO Xce allele. Now we want to compare to B6 mutant allele (XceB) and see which is stronger
#only look at females
Xce_allele <- Xce_allele[rownames(pheno.F),]
Xce_allele <- as.data.frame( Xce_allele)
names(Xce_allele) <- "Xce_DO"
Xce_allele$Xce_Mut <- "Xce_b"
Xce_allele$Xce_winner <- 0

sample_names <- rownames(Xce_allele)
for (sample in sample_names){

	if ((Xce_allele[sample, "Xce_DO"] == "Xce_a") || (Xce_allele[sample, "Xce_DO"] == "Xce_b") ||
		(Xce_allele[sample, "Xce_DO"] == "Xce_e")){

		Xce_allele[sample, "Xce_winner"] <- "Xce_b"
	}
	else{
		Xce_allele[sample, "Xce_winner"] <- "Xce_c"
	}
}

#lable origin either: DO, Mut, Both
Xce_allele$Winner_origin <- 0

for (sample in sample_names){

	if ((Xce_allele[sample, "Xce_winner"] == Xce_allele[sample, "Xce_DO"]) & (Xce_allele[sample, "Xce_winner"] == Xce_allele[sample, "Xce_Mut"])){
		Xce_allele[sample, "Winner_origin"] <- "Both"
	}
	else if (Xce_allele[sample, "Xce_winner"] == Xce_allele[sample, "Xce_DO"]){
		Xce_allele[sample, "Winner_origin"] <- "DO"
	}
	else if (Xce_allele[sample, "Xce_winner"] == Xce_allele[sample, "Xce_Mut"]){
		Xce_allele[sample, "Winner_origin"] <- "Mut"
	}
}

save(Xce_allele, file = paste0(dir,"Data/consolidated/Phenotype/Xce_allele.Rdata"))

####################################################################################################################################################################
####################################################################################################################################################################
##Xist
#	In this section we will be using the output from EMASE in determining which founder contributed to the Xist expression
#	Luckily Xist is already annotated and the gene ID name is known so it should be a simple look up wiht a little bit of cleaning
#	Xist Ensembl gene ID is : ENSMUSG00000086503

#The following is done on the Cadillac linux commandline
cd /hpcdata/ytakemon/Col4a5xDO/civet_run
grep -H ENSMUSG00000086503 ./*/gbrs.quantified.diploid.genes.tpm > /hpcdata/ytakemon/Col4a5xDO/Phenotype/Xist_diploid_exp.txt

#Back to R
Xist <- read.delim(file = "/hpcdata/ytakemon/Col4a5xDO/Phenotype/Xist_diploid_exp.txt", sep = "\t", header = TRUE)
rownames(Xist) <- Xist$MouseID
Xist$MouseID <- NULL
#The expression varies greatly so calculte the percentage form the total for each founder
Xist$total <- rowSums(Xist[,c(1:8)])
Xist[,1:9] <- Xist[,1:9] / Xist[,"total"]

#find founder with the max % of XIST expression
sample_names <- rownames(Xist)
Xist_max_founder <- array(0, c(192, 2), dimnames  = list( sample_names, c("max_founders", "max_perc"))) #192 samples, 2col(allele, maxprob),
for (sample in sample_names){

	df <- Xist[sample, c(1:8)]
	max.founder <- colnames(Xist[,c(1:8)])[apply(df,1,which.max)] #2 for column, #1 for row
	max.prob <- max(df)

	Xist_max_founder[sample, "max_founders" ] <- max.founder
	Xist_max_founder[sample,	"max_perc" ] <- max.prob
}
Xist_max_founder <- as.data.frame(Xist_max_founder)
Xist_max_founder$dip_called <- Xist$genotype

save(Xist_max_founder, file = "/hpcdata/ytakemon/Col4a5xDO/Phenotype/Max_Xist_founder.Rdata")
####################################################################################################################################################################
####################################################################################################################################################################
## Compare Xce with Xist expression
load("/hpcdata/ytakemon/Col4a5xDO/Phenotype/Max_Xist_founder.Rdata")
load("./Phenotype/Xce_allele.Rdata")

#subset Xist to only females
Xist_max_founder <- Xist_max_founder[rownames(Xce_allele),]

##	HYPOTHESIS
# if dominant Xce is coming from the Mut (B6), then the genotype silenced (Xist) comes from the DO (A-H).
# if dominant Xce is coming from the DO , then the genotype silenced must be B.
# if dominnat Xce comes from both , them the genotype silenced is 50/50 cannot be determined.

#test
Xce_allele$Xist_maxfounder <- Xist_max_founder$max_founders
Xce_allele$Hypothesis_test <- 0

sample_names <- rownames(Xce_allele)
for (sample in sample_names){

	Test <- Xce_allele[sample, "Winner_origin"]
	Xist <- Xce_allele[sample, "Xist_maxfounder"]

	if (Test == "Mut"){
		Xce_allele[sample, "Hypothesis_test"] <- "TRUE"
	}
	else if (Test == "DO"){
		if (Xist == "B"){
			Xce_allele[sample, "Hypothesis_test"] <- "TRUE"
		}
		else{
			Xce_allele[sample, "Hypothesis_test"] <- "FALSE"
		}
	}
	else if (Test == "Both"){
		Xce_allele[sample, "Hypothesis_test"] <- "TRUE"
	}

}

#Hypothesis holds true mostly except for the following samples
# Sample 	Dom.Xce 	Xist
#1415.0201	DO 			F
#1415.0813	DO 			A
#1415.0934	DO 			H
#1415.0977	DO 			F

#Xce silencing probability is only looked at in Inbred mice, so when the cis Xist initiation site is of differnt origin, it isn't as effective?
#also this is just a probability.


#########################################ANOVA TEST
##ANOVA test between Xce allele and phenotype
#	Test to see if phenotypes such as ACR and GFR are assocaited with the orgigin of Xce allele
#load and organize data
load(paste0(dir,"Data/consolidated/Phenotype/Xce_allele.Rdata"))
pheno <- read.delim(paste0(dir,"Data/consolidated/Phenotype/1415_master_pheno.txt"), sep = "\t", header = TRUE)
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(Xce_allele),] #subset pheno to match only females

pheno$Winner_origin <- Xce_allele$Winner_origin
pheno$Xce_DO <- Xce_allele$Xce_DO
pheno$Xce_winner <- Xce_allele$Xce_winner

#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA

pheno$C2_log <- log(pheno$C2)

pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR10WK_log <- log(pheno$ACR10WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)

pheno$Alb6WK_log <- log(pheno$Alb6WK)
pheno$Alb10WK_log <- log(pheno$Alb10WK)
pheno$Alb15WK_log <- log(pheno$Alb15WK)

pheno$Creat6WK_log <- log(pheno$Creat6WK)
pheno$Creat10WK_log <- log(pheno$Creat10WK)
pheno$Creat15WK_log <- log(pheno$Creat15WK)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs

#remove unncessary columns
pheno <- pheno[,c("MouseID", "C2_log", "Alb6WK_log","Creat6WK_log","Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log", "Winner_origin", "Xce_DO", "Xce_winner")]

####
table(pheno$Xce_winner)
table(pheno$Winner_origin)

###### ANOVA test between Xce_winner and phenotypes
#Actually Xce_winner only has 2, so T-test.... wiat acutally ANOVA because of covariates
#first plot them
library(ggplot2)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_winner_GFR.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_winner, y = C2_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce alleles") +
	ylab( "Log-transformed C2 GFR") +
	ggtitle("Box plot of Winning Xce allele v.s. log(GFR)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_winner_Alb6.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_winner, y = Alb6WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce alleles") +
	ylab( "Log-transformed Alb6wk") +
	ggtitle("Box plot of Winning Xce allele v.s. log(Albumin 6WK")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_winner_Alb10.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_winner, y = Alb10WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce alleles") +
	ylab( "Log-transformed Alb10wk") +
	ggtitle("Box plot of Winning Xce allele v.s. log(Albumin 10WK")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_winner_Alb15.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_winner, y = Alb15WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce alleles") +
	ylab( "Log-transformed Alb15wk") +
	ggtitle("Box plot of Winning Xce allele v.s. log(Albumin 15WK")
dev.off()

#ANOVA: Xce_winner vs. Phenotype
##one-way ANOVA Test format:
#	fit <- aov(y ~ A + x, data=mydataframe)
#	y : numeric value to test
#	A : factors/ categories
#	x : if there is a covariate use x
##note:
#	aov() and lm() gives different outputs, however if the object is then run with anova() the resutls are the same.

##GFR ANOVA Xce_winner vs. Phenotype
G.pheno <- pheno[complete.cases(pheno$C2_log),]
GFR_lm <- lm(C2_log ~ Xce_winner, data = G.pheno)
GFR_anova <- aov(C2_log ~ Xce_winner, data = G.pheno)
anova(GFR_lm)
knitr::kable(anova(GFR_anova))
##ANOVA result
#no sign difference
#ANOVA test assumptions for heteroscadacity
GFR_lm_resid <- data.frame(Fitted = fitted(GFR_lm), Residuals = resid(GFR_lm), Xce = G.pheno$Xce_winner)
#plot residuals
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_winner_GFR.png", width = 1500, height = 1000, res = 100)
ggplot(GFR_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_winner_GFR.png", width = 1500, height = 1000, res = 100)
plot(GFR_anova)
dev.off()
# Residuals look normally distributed

##6wk Albumin ANOVA Xce_winner vs. Phenotype
#subset data
A6.pheno <- pheno[complete.cases(pheno$Alb6WK_log),]
A6.pheno <- A6.pheno[complete.cases(A6.pheno$Creat6WK_log),]
A6_lm <- lm(Alb6WK_log ~ Xce_winner, data = A6.pheno)
A6_anova <- aov(Alb6WK_log ~ Xce_winner, data = A6.pheno)
A6_covar_lm <- lm(Alb6WK_log ~ Xce_winner + Creat6WK_log, data = A6.pheno)
A6_covar_anova <- aov(Alb6WK_log ~ Xce_winner + Creat6WK_log, data = A6.pheno)
anova(A6_lm)
knitr::kable(anova(A6_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A6_lm, A6_covar_lm) #different
## ANOVA results
#there is a significant difference between groups when creatinine is considered a covariate
#plot residuals
A6_covar_lm_resid <- data.frame(Fitted = fitted(A6_covar_lm), Residuals = resid(A6_covar_lm), Xce = A6.pheno$Xce_winner)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_winner_Alb6.png", width = 1500, height = 1000, res = 100)
ggplot(A6_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_winner_Alb6.png", width = 1500, height = 1000, res = 100)
plot(A6_covar_anova)
dev.off()

##10wk Albumin ANOVA Xce_winner vs. Phenotype
#subset data
A10.pheno <- pheno[complete.cases(pheno$Alb10WK_log),]
A10.pheno <- A10.pheno[complete.cases(A10.pheno$Creat10WK_log),]
A10_lm <- lm(Alb10WK_log ~ Xce_winner, data = A10.pheno)
A10_anova <- aov(Alb10WK_log ~ Xce_winner, data = A10.pheno)
A10_covar_lm <- lm(Alb10WK_log ~ Xce_winner + Creat10WK_log, data = A10.pheno)
A10_covar_anova <- aov(Alb10WK_log ~ Xce_winner + Creat10WK_log, data = A10.pheno)
anova(A10_lm)
knitr::kable(anova(A10_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A10_lm, A10_covar_lm) #not different
## ANOVA results
#there is a significant difference between groups when creatinine is considered a covariate
#plot residuals
A10_covar_lm_resid <- data.frame(Fitted = fitted(A10_covar_lm), Residuals = resid(A10_covar_lm), Xce = A10.pheno$Xce_winner)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_winner_Alb10.png", width = 1500, height = 1000, res = 100)
ggplot(A10_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_winner_Alb10.png", width = 1500, height = 1000, res = 100)
plot(A10_covar_anova)
dev.off()

##15wk Albumin ANOVA Xce_winner vs. Phenotype
#subset data
A15.pheno <- pheno[complete.cases(pheno$Alb15WK_log),]
A15.pheno <- A15.pheno[complete.cases(A15.pheno$Creat15WK_log),]
A15_lm <- lm(Alb15WK_log ~ Xce_winner, data = A15.pheno)
A15_anova <- aov(Alb15WK_log ~ Xce_winner, data = A15.pheno)
A15_covar_lm <- lm(Alb15WK_log ~ Xce_winner + Creat15WK_log, data = A15.pheno)
A15_covar_anova <- aov(Alb15WK_log ~ Xce_winner + Creat15WK_log, data = A15.pheno)
anova(A15_lm)
knitr::kable(anova(A15_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A15_lm, A15_covar_lm) #different
## ANOVA results
#there is no significant difference
#plot residuals
A15_covar_lm_resid <- data.frame(Fitted = fitted(A15_covar_lm), Residuals = resid(A15_covar_lm), Xce = A15.pheno$Xce_winner)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_winner_Alb15.png", width = 1500, height = 1000, res = 100)
ggplot(A15_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_winner_Alb15.png", width = 1500, height = 1000, res = 100)
plot(A15_covar_anova)
dev.off()


#ANOVA: Winner_origin vs. Phenotype
##one-way ANOVA Test format:
#	fit <- aov(y ~ A + x, data=mydataframe)
#	y : numeric value to test
#	A : factors/ categories
#	x : if there is a covariate use x
##note:
#	aov() and lm() gives different outputs, however if the object is then run with anova() the resutls are the same.

##Make boxplot for Winner_origin vs. Phenotype
library(ggplot2)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Winner_origin_GFR.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Winner_origin, y = C2_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce allele origin") +
	ylab( "Log-transformed C2 GFR") +
	ggtitle("Box plot of Winning allele origin v.s. log(GFR)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Winner_origin_Alb6.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Winner_origin, y = Alb6WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce allele origin") +
	ylab( "Log-transformed Alb6wk")+
	ggtitle("Box plot of Winning allele origin v.s. log(Albumin 6WK)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Winner_origin_Alb10.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Winner_origin, y = Alb10WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce allele origin") +
	ylab( "Log-transformed Alb10wk")+
	ggtitle("Box plot of Winning allele origin v.s. log(Albumin 10WK)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Winner_origin_Alb15.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Winner_origin, y = Alb15WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "Xce allele origin") +
	ylab( "Log-transformed Alb15wk")+
	ggtitle("Box plot of Winning allele origin v.s. log(Albumin 15WK)")
dev.off()

##GFR ANOVA: Winner_origin vs. Phenotype
G.pheno <- pheno[complete.cases(pheno$C2_log),]
GFR_lm <- lm(C2_log ~ Winner_origin, data = G.pheno)
GFR_anova <- aov(C2_log ~ Winner_origin, data = G.pheno)
anova(GFR_lm)
knitr::kable((anova(GFR_anova)))
##ANOVA result
#no sign difference
#ANOVA test assumptions for heteroscadacity
GFR_lm_resid <- data.frame(Fitted = fitted(GFR_lm), Residuals = resid(GFR_lm), Xce = G.pheno$Winner_origin)
#plot residuals
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Winner_origin_GFR.png", width = 1500, height = 1000, res = 100)
ggplot(GFR_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Winner_origin_GFR.png", width = 1500, height = 1000, res = 100)
plot(GFR_anova)
dev.off()
# Residuals look normally distributed
##	Tukey's HSD test
posthoc_G <- TukeyHSD( x = GFR_anova, "Winner_origin", conf.level = 0.95)
kable(posthoc_G$Winner_origin)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Winner_origin_GFR.png", width = 1500, height = 1000, res = 100)
plot(posthoc_G)
dev.off()

##6wk Albumin ANOVA : Winner_origin vs. Phenotype
#subset data
A6.pheno <- pheno[complete.cases(pheno$Alb6WK_log),]
A6.pheno <- A6.pheno[complete.cases(A6.pheno$Creat6WK_log),]
A6_lm <- lm(Alb6WK_log ~ Winner_origin, data = A6.pheno)
A6_anova <- aov(Alb6WK_log ~ Winner_origin, data = A6.pheno)
A6_covar_lm <- lm(Alb6WK_log ~ Winner_origin + Creat6WK_log, data = A6.pheno)
A6_covar_anova <- aov(Alb6WK_log ~ Winner_origin + Creat6WK_log, data = A6.pheno)
anova(A6_lm)
knitr::kable(anova(A6_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A6_lm, A6_covar_lm) #different
## ANOVA results
#there is a significant difference between groups when creatinine is considered a covariate
#plot residuals
A6_covar_lm_resid <- data.frame(Fitted = fitted(A6_covar_lm), Residuals = resid(A6_covar_lm), Xce = A6.pheno$Winner_origin)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Winner_origin_Alb6.png", width = 1500, height = 1000, res = 100)
ggplot(A6_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Winner_origin_Alb6.png", width = 1500, height = 1000, res = 100)
plot(A6_covar_anova)
dev.off()
##	Tukey's HSD test
posthoc_A6 <- TukeyHSD( x = A6_covar_anova, "Winner_origin", conf.level = 0.95)
kable(posthoc_A6$Winner_origin)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Winner_origin_Alb6.png", width = 1500, height = 1000, res = 100)
plot(posthoc_A6)
dev.off()

##10wk Albumin ANOVA : Winner_origin vs. Phenotype
#subset data
A10.pheno <- pheno[complete.cases(pheno$Alb10WK_log),]
A10.pheno <- A10.pheno[complete.cases(A10.pheno$Creat10WK_log),]
A10_lm <- lm(Alb10WK_log ~ Winner_origin, data = A10.pheno)
A10_anova <- aov(Alb10WK_log ~ Winner_origin, data = A10.pheno)
A10_covar_lm <- lm(Alb10WK_log ~ Winner_origin + Creat10WK_log, data = A10.pheno)
A10_covar_anova <- aov(Alb10WK_log ~ Winner_origin + Creat10WK_log, data = A10.pheno)
anova(A10_lm)
knitr::kable((anova(A10_covar_lm)))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A10_lm, A10_covar_lm) #not different
## ANOVA results
#there is a significant difference between groups when creatinine is considered a covariate
#plot residuals
A10_covar_lm_resid <- data.frame(Fitted = fitted(A10_covar_lm), Residuals = resid(A10_covar_lm), Xce = A10.pheno$Winner_origin)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Winner_origin_Alb10.png", width = 1500, height = 1000, res = 100)
ggplot(A10_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Winner_origin_Alb10.png", width = 1500, height = 1000, res = 100)
plot(A10_covar_anova)
dev.off()
##	Tukey's HSD test
posthoc_A10 <- TukeyHSD( x = A10_covar_anova, "Winner_origin", conf.level = 0.95)
kable(posthoc_A10$Winner_origin)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Winner_origin_Alb10.png", width = 1500, height = 1000, res = 100)
plot(posthoc_A10)
dev.off()

##15wk Albumin ANOVA : Winner_origin vs. Phenotype
#subset data
A15.pheno <- pheno[complete.cases(pheno$Alb15WK_log),]
A15.pheno <- A15.pheno[complete.cases(A15.pheno$Creat15WK_log),]
A15_lm <- lm(Alb15WK_log ~ Winner_origin, data = A15.pheno)
A15_anova <- aov(Alb15WK_log ~ Winner_origin, data = A15.pheno)
A15_covar_lm <- lm(Alb15WK_log ~ Winner_origin + Creat15WK_log, data = A15.pheno)
A15_covar_anova <- aov(Alb15WK_log ~ Winner_origin + Creat15WK_log, data = A15.pheno)
anova(A15_lm)
knitr::kable((anova(A15_covar_lm)))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A15_lm, A15_covar_lm) #not different
## ANOVA results
#there is no difference
#plot residuals
A15_covar_lm_resid <- data.frame(Fitted = fitted(A15_covar_lm), Residuals = resid(A15_covar_lm), Xce = A15.pheno$Winner_origin)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Winner_origin_Alb15.png", width = 1500, height = 1000, res = 100)
ggplot(A15_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Winner_origin_Alb15.png", width = 1500, height = 1000, res = 100)
plot(A15_covar_anova)
dev.off()
##	Tukey's HSD test
posthoc_A15 <- TukeyHSD( x = A15_covar_anova, "Winner_origin", conf.level = 0.95)
kable(posthoc_A15$Winner_origin)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Winner_origin_Alb15.png", width = 1500, height = 1000, res = 100)
plot(posthoc_A15)
dev.off()


#ANOVA: Xce_DO vs. Phenotype
##one-way ANOVA Test format:
#	fit <- aov(y ~ A + x, data=mydataframe)
#	y : numeric value to test
#	A : factors/ categories
#	x : if there is a covariate use x
##note:
#	aov() and lm() gives different outputs, however if the object is then run with anova() the resutls are the same.

##Make avg with std.error for Xce alleles vs ACR 6, 10, 15
data <- pheno %>%
  mutate(sample = rownames(pheno)) %>%
  select(sample, Xce_DO, ACR6WK_log, ACR10WK_log, ACR15WK_log)

data6 <- data[,c(1,2)]
data10 <- data[,c(1,3)]
data15 <- data[,c(1,4)]

data6log <- data[,c(1,5)]
data10log <- data[,c(1,6)]
data15log <- data[,c(1,7)]

ggdata6 <- summarySE(data6, measurevar = "ACR6WK", groupvars = "Xce_DO", na.rm = TRUE)
ggdata10 <- summarySE(data10, measurevar = "ACR10WK", groupvars = "Xce_DO", na.rm = TRUE)
ggdata15 <- summarySE(data15, measurevar = "ACR15WK", groupvars = "Xce_DO", na.rm = TRUE)
ggdata6log <- summarySE(data6log, measurevar = "ACR6WK_log", groupvars = "Xce_DO", na.rm = TRUE)
ggdata10log <- summarySE(data10log, measurevar = "ACR10WK_log", groupvars = "Xce_DO", na.rm = TRUE)
ggdata15log <- summarySE(data15log, measurevar = "ACR15WK_log", groupvars = "Xce_DO", na.rm = TRUE)


ggdata <- data %>%
  pivot_longer(-c(sample,Xce_DO), names_to = "ACR_month", values_to = "Value")
p_load(ggsci)

png(paste0(dir,"Results/Xce_allele_ACR.png"), height = 4, width = 6, units = "in", compression = "none", res = )
ggdata %>%
  mutate(Month = case_when(
            .$ACR_month == "ACR6WK_log" ~ "6",
            .$ACR_month == "ACR10WK_log" ~ "10",
            .$ACR_month == "ACR15WK_log" ~ "15",
            TRUE ~ "NA")) %>%
  mutate(Month = fct_relevel(Month, "6","10","15")) %>%
ggplot(., aes(x = Month, y = Value, fill = Xce_DO))+
  geom_boxplot() +
  scale_fill_aaas() +
  labs(x = "Weeks of age",
       y = "log(ACR)",
       fill = "Xce allele") +
  theme_bw()+
  ylim(1,10)
dev.off()

test <- ggdata %>%
  mutate(Month = case_when(
            .$ACR_month == "ACR6WK_log" ~ "6",
            .$ACR_month == "ACR10WK_log" ~ "10",
            .$ACR_month == "ACR15WK_log" ~ "15",
            TRUE ~ "NA")) %>%
  mutate(Month = fct_relevel(Month, "6","10","15"))

# repeated measures ANOVA in R
test %>%
  filter(complete.cases(Value)) %>%
  aov(Value ~ Xce_DO + Error(sample/Month), data = .) %>% tidy()
# A tibble: 3 x 7
#  stratum      term         df sumsq meansq statistic  p.value
#  <chr>        <chr>     <dbl> <dbl>  <dbl>     <dbl>    <dbl>
#1 sample       Xce_DO        3  36.8  12.3       5.18  0.00236
#2 sample       Residuals    92 218.    2.37     NA    NA
#3 sample:Month Residuals   152 241.    1.59     NA    NA


# post hoc
test %>%
  filter(Month == "10") %>%
  filter(complete.cases(Value)) %>%
  aov(Value ~ Xce_DO, + Error().) %>% tidy()

# post hoc
test %>%
  filter(Month == "10") %>%
  filter(complete.cases(Value)) %>%
  aov(Value ~ Xce_DO, .) %>%
  TukeyHSD(.,"Xce_DO")

source("~/GitHub/1415-Col4a5xDO-Project/RHelper_func/STDEV_SE_NORM_functions.R")
ggdata <- summarySE(ggdata, measurevar = "Value", groupvars = c("Xce_DO", "ACR_month"), na.rm = TRUE)

allggplot <- ggplot(ggdata, aes(x = ACR_month, y = Value, colour = Xce_DO)) +
	geom_line(aes( group = Xce_DO)) +
	geom_errorbar(aes(ymin = Value - se, ymax = Value + se), width = 0.2) +
	geom_point() +
	labs(title = "Log-transformed ACR by DO Xce Alleles", x = "Months", y = "Log-transformed ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))


ggplot6 <- ggplot(ggdata6, aes(x = Xce_DO, y = ACR6WK, colour = Xce_DO)) +
	geom_errorbar(aes( ymin = ACR6WK - se, ymax = ACR6WK + se), width = 0.1) +
	geom_point() +
	labs( title = "6-week ACR by DO Xce Alleles", x = "DO Xce alleles", y = "6-week ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

ggplot10 <- ggplot(ggdata10, aes(x = Xce_DO, y = ACR10WK, colour = Xce_DO)) +
	geom_errorbar(aes( ymin = ACR10WK - se, ymax = ACR10WK + se), width = 0.1) +
	geom_point() +
	labs( title = "10-week ACR by DO Xce Alleles", x = "DO Xce alleles", y = "10-week ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

ggplot15 <- ggplot(ggdata15, aes(x = Xce_DO, y = ACR15WK, colour = Xce_DO)) +
	geom_errorbar(aes( ymin = ACR15WK - se, ymax = ACR15WK + se), width = 0.1) +
	geom_point() +
	labs( title = "15-week ACR by DO Xce Alleles", x = "DO Xce alleles", y = "15-week ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

ggplot6log <- ggplot(ggdata6log, aes(x = Xce_DO, y = ACR6WK_log, colour = Xce_DO)) +
	geom_errorbar(aes( ymin = ACR6WK_log - se, ymax = ACR6WK_log + se), width = 0.1) +
	geom_point() +
	scale_y_continuous(limits = c(4.0,8.0)) +
	labs( title = "Log-transformed 6-week ACR by DO Xce Alleles", x = "DO Xce alleles", y = "6-week ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

ggplot10log <- ggplot(ggdata10log, aes(x = Xce_DO, y = ACR10WK_log, colour = Xce_DO)) +
	geom_errorbar(aes( ymin = ACR10WK_log - se, ymax = ACR10WK_log + se), width = 0.1) +
	geom_point() +
	scale_y_continuous(limits = c(4.0,8.0)) +
	labs( title = "Log-transformed 10-week ACR by DO Xce Alleles", x = "DO Xce alleles", y = "10-week ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

ggplot15log <- ggplot(ggdata15log, aes(x = Xce_DO, y = ACR15WK_log, colour = Xce_DO)) +
	geom_errorbar(aes( ymin = ACR15WK_log - se, ymax = ACR15WK_log + se), width = 0.1) +
	geom_point() +
	scale_y_continuous(limits = c(4.0,8.0)) +
	labs( title = "Log-transformed 15-week ACR by DO Xce Alleles", x = "DO Xce alleles", y = "15-week ACR") +
	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_All_logACR.pdf", width = 10.0, height = 7.5)
allggplot
dev.off()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_ACR6WK.pdf", width = 10.0, height = 7.5)
ggplot6
dev.off()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_ACR10WK.pdf", width = 10.0, height = 7.5)
ggplot10
dev.off()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_ACR15WK.pdf", width = 10.0, height = 7.5)
ggplot15
dev.off()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_ACR6WK_log.pdf", width = 10.0, height = 7.5)
ggplot6log
dev.off()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_ACR10WK_log.pdf", width = 10.0, height = 7.5)
ggplot10log
dev.off()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_stderror_Xce_DO_ACR15WK_log.pdf", width = 10.0, height = 7.5)
ggplot15log
dev.off()

#view statistics
ggdata6 <- data6[complete.cases(data6$ACR6WK),]
ggdata10 <- data10[complete.cases(data10$ACR10WK),]
ggdata15 <- data15[complete.cases(data15$ACR15WK),]

ggdata6log <- data6log[complete.cases(data6log$ACR6WK),]
ggdata10log <- data10log[complete.cases(data10log$ACR10WK),]
ggdata15log <- data15log[complete.cases(data15log$ACR15WK),]

lm_6 <- aov(ACR6WK ~ Xce_DO, data = ggdata6)
kable(anova(lm_6))
lm_10 <- aov(ACR10WK ~ Xce_DO, data = ggdata10)
kable(anova(lm_10))
lm_15 <- aov(ACR15WK ~ Xce_DO, data = ggdata15)
kable(anova(lm_15))

lm_6_log <- aov(ACR6WK_log ~ Xce_DO, data = ggdata6log)
kable(anova(lm_6_log))
lm_10_log <- aov(ACR10WK_log ~ Xce_DO, data = ggdata10log)
kable(anova(lm_10_log))
lm_15_log <- aov(ACR15WK_log ~ Xce_DO, data = ggdata15log)
kable(anova(lm_15_log))

posthoc_6log <- TukeyHSD( x = lm_6_log, "Xce_DO", conf.level = 0.95)
kable(posthoc_6log$Xce_DO)
posthoc_10log <- TukeyHSD( x = lm_10_log, "Xce_DO", conf.level = 0.95)
kable(posthoc_10log$Xce_DO)
posthoc_15log <- TukeyHSD( x = lm_15_log, "Xce_DO", conf.level = 0.95)
kable(posthoc_15log$Xce_DO)


##Make boxplot for Xce_DO vs. Phenotype
library(ggplot2)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_DO_GFR.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_DO, y = C2_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "DO Xce allele") +
	ylab( "Log-transformed C2 GFR") +
	ggtitle("Box plot of DO Xce allele v.s. log(GFR)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_DO_Alb6.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_DO, y = Alb6WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "DO Xce allele") +
	ylab( "Log-transformed Alb6wk") +
	ggtitle("Box plot of DO Xce allele v.s. log(Albumin 6WK)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_DO_Alb10.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_DO, y = Alb10WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "DO Xce allele") +
	ylab( "Log-transformed Alb10wk") +
	ggtitle("Box plot of DO Xce allele v.s. log(Albumin 10WK)")
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Boxplot_Xce_DO_Alb15.png", width = 1500, height = 1000, res = 100)
ggplot( pheno, aes(x  = Xce_DO, y = Alb15WK_log)) +
	geom_boxplot( fill = "grey80", colour = "blue") +
	scale_x_discrete() + theme(text = element_text(size =25)) +
	xlab( "DO Xce allele") +
	ylab( "Log-transformed Alb15wk") +
	ggtitle("Box plot of DO Xce allele v.s. log(Albumin 15WK)")
dev.off()

##GFR ANOVA: Xce_DO vs. Phenotype
G.pheno <- pheno[complete.cases(pheno$C2_log),]
GFR_lm <- lm(C2_log ~ Xce_DO, data = G.pheno)
GFR_anova <- aov(C2_log ~ Xce_DO, data = G.pheno)
anova(GFR_lm)
kable(anova(GFR_anova))
##ANOVA result
#no sign difference
#ANOVA test assumptions for heteroscadacity
GFR_lm_resid <- data.frame(Fitted = fitted(GFR_lm), Residuals = resid(GFR_lm), Xce = G.pheno$Xce_DO)
#plot residuals
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_DO_GFR.png", width = 1500, height = 1000, res = 100)
ggplot(GFR_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_DO_GFR.png", width = 1500, height = 1000, res = 100)
plot(GFR_anova)
dev.off()
# Residuals look normally distributed
##	Tukey's HSD test
posthoc_G <- TukeyHSD( x = GFR_anova, "Xce_DO", conf.level = 0.95)
kable(posthoc_G$Xce_DO)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Xce_DO_GFR.png", width = 1500, height = 1000, res = 100)
plot(posthoc)
dev.off()

##6wk Albumin ANOVA : Xce_DO vs. Phenotype
#subset data
A6.pheno <- pheno[complete.cases(pheno$Alb6WK_log),]
A6.pheno <- A6.pheno[complete.cases(A6.pheno$Creat6WK_log),]
A6_lm <- lm(Alb6WK_log ~ Xce_DO, data = A6.pheno)
A6_anova <- aov(Alb6WK_log ~ Xce_DO, data = A6.pheno)
A6_covar_lm <- lm(Alb6WK_log ~ Xce_DO + Creat6WK_log, data = A6.pheno)
A6_covar_anova <- aov(Alb6WK_log ~ Xce_DO + Creat6WK_log, data = A6.pheno)
anova(A6_lm)
kable(anova(A6_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A6_lm, A6_covar_lm) #different
## ANOVA results
#there is a significant difference between groups when creatinine is considered a covariate
#plot residuals
A6_covar_lm_resid <- data.frame(Fitted = fitted(A6_covar_lm), Residuals = resid(A6_covar_lm), Xce = A6.pheno$Xce_DO)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_DO_Alb6.png", width = 1500, height = 1000, res = 100)
ggplot(A6_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_DO_Alb6.png", width = 1500, height = 1000, res = 100)
plot(A6_covar_anova)
dev.off()
##	Tukey's HSD test
posthoc_A6 <- TukeyHSD( x = A6_covar_anova, "Xce_DO", conf.level = 0.95)
kable(posthoc_A6$Xce_DO)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Xce_DO_Alb6.png", width = 1500, height = 1000, res = 100)
plot(posthoc_A6)
dev.off()

##10wk Albumin ANOVA : Xce_DO vs. Phenotype
#subset data
A10.pheno <- pheno[complete.cases(pheno$Alb10WK_log),]
A10.pheno <- A10.pheno[complete.cases(A10.pheno$Creat10WK_log),]
A10_lm <- lm(Alb10WK_log ~ Xce_DO, data = A10.pheno)
A10_anova <- aov(Alb10WK_log ~ Xce_DO, data = A10.pheno)
A10_covar_lm <- lm(Alb10WK_log ~ Xce_DO + Creat10WK_log, data = A10.pheno)
A10_covar_anova <- aov(Alb10WK_log ~ Xce_DO + Creat10WK_log, data = A10.pheno)
anova(A10_lm)
kable(anova(A10_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A10_lm, A10_covar_lm) #different
## ANOVA results
#there is a significant difference between groups when creatinine is considered a covariate
#plot residuals
A10_covar_lm_resid <- data.frame(Fitted = fitted(A10_covar_lm), Residuals = resid(A10_covar_lm), Xce = A10.pheno$Xce_DO)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_DO_Alb10.png", width = 1500, height = 1000, res = 100)
ggplot(A10_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_DO_Alb10.png", width = 1500, height = 1000, res = 100)
plot(A10_covar_anova)
dev.off()
##	Tukey's HSD test
posthoc_A10 <- TukeyHSD( x = A10_covar_anova, "Xce_DO", conf.level = 0.95)
kable(posthoc_A10$Xce_DO)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Xce_DO_Alb10.png", width = 1500, height = 1000, res = 100)
plot(posthoc_A10)
dev.off()

##15wk Albumin ANOVA : Xce_DO vs. Phenotype
#subset data
A15.pheno <- pheno[complete.cases(pheno$Alb15WK_log),]
A15.pheno <- A15.pheno[complete.cases(A15.pheno$Creat15WK_log),]
A15_lm <- lm(Alb15WK_log ~ Xce_DO, data = A15.pheno)
A15_anova <- aov(Alb15WK_log ~ Xce_DO, data = A15.pheno)
A15_covar_lm <- lm(Alb15WK_log ~ Xce_DO + Creat15WK_log, data = A15.pheno)
A15_covar_anova <- aov(Alb15WK_log ~ Xce_DO + Creat15WK_log, data = A15.pheno)
anova(A15_lm)
kable(anova(A15_covar_lm))
#compare both outcomes to see if creatinine covariate created significant difference in results
anova(A15_lm, A15_covar_lm) #not different
## ANOVA results
#there is no difference
#plot residuals
A15_covar_lm_resid <- data.frame(Fitted = fitted(A15_covar_lm), Residuals = resid(A15_covar_lm), Xce = A15.pheno$Xce_DO)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual_Xce_DO_Alb15.png", width = 1500, height = 1000, res = 100)
ggplot(A15_covar_lm_resid, aes( Fitted, Residuals, colour = Xce)) + theme(text = element_text(size =25)) + geom_point()
dev.off()
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Residual2_Xce_DO_Alb15.png", width = 1500, height = 1000, res = 100)
plot(A15_covar_anova)
dev.off()
##	Tukey's HSD test
posthoc_A15 <- TukeyHSD( x = A15_covar_anova, "Xce_DO", conf.level = 0.95)
kable(posthoc_A15$Xce_DO)
#plot Tukey's HSD test
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Tuckey_Xce_Xce_DO_Alb15.png", width = 1500, height = 1000, res = 100)
plot(posthoc_A15)
dev.off()
