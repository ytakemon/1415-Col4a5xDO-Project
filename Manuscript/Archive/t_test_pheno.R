library(reshape2)
library(ggplot2)
library(grid)

setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
#keep only interested columns
pheno <- pheno[,c("MouseID", "Sex", "C2", "ACR6WK", "ACR10WK", "ACR15WK")]

pheno_F <- pheno[pheno$Sex == "F",]
pheno_M <- pheno[pheno$Sex == "M",]

#remove NAs from C2 data
C2_F <- pheno_F$C2[complete.cases(pheno_F$C2)]
C2_M <- pheno_M$C2[complete.cases(pheno_M$C2)]
ttest_C2 <- t.test(C2_F, C2_M)
# 	Welch Two Sample t-test
#
# data:  pheno_F$C2 and pheno_M$C2
# t = 2.1992, df = 85.347, p-value = 0.03057
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#    6.137505 121.753695
# sample estimates:
# mean of x mean of y
#  286.8099  222.8643
#
# standard error of mean
sem_C2_F<-sd(C2_F)/sqrt(length(C2_F)) # 12.32988
sem_C2_M<-sd(C2_M)/sqrt(length(C2_M)) # 26.3326

# remove NAs from the ACR6WK data
ACR6WK_F <- pheno_F$ACR6WK[complete.cases(pheno_F$ACR6WK)]
ACR6WK_M <- pheno_M$ACR6WK[complete.cases(pheno_M$ACR6WK)]
ttest_ACR6WK <- t.test(ACR6WK_F, ACR6WK_M)
#	Welch Two Sample t-test
#
# data:  ACR6WK_F and ACR6WK_M
# t = -3.7029, df = 116.92, p-value = 0.0003268
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -860.2846 -260.7203
# sample estimates:
# mean of x mean of y
#  427.7154  988.2179
#
# standard error of mean
sem_ACR6WK_F<-sd(ACR6WK_F)/sqrt(length(ACR6WK_F)) # 62.6336
sem_ACR6WK_M<-sd(ACR6WK_M)/sqrt(length(ACR6WK_M)) # 137.8037

# remove NAs from the ACR10WK data
ACR10WK_F <- pheno_F$ACR10WK[complete.cases(pheno_F$ACR10WK)]
ACR10WK_M <- pheno_M$ACR10WK[complete.cases(pheno_M$ACR10WK)]
ttest_ACR10WK <- t.test(ACR10WK_F, ACR10WK_M)
#	Welch Two Sample t-test
#
# data:  ACR10WK_F and ACR10WK_M
# t = -6.8251, df = 112.59, p-value = 4.66e-10
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -2617.652 -1439.811
# sample estimates:
# mean of x mean of y
#  1581.270  3610.001
#
# standard error of mean
sem_ACR10WK_F<-sd(ACR10WK_F)/sqrt(length(ACR10WK_F)) # 133.4909
sem_ACR10WK_M<-sd(ACR10WK_M)/sqrt(length(ACR10WK_M)) # 265.5846

#remove NAs from ACR15WK data
ACR15WK_F <- pheno_F$ACR15WK[complete.cases(pheno_F$ACR15WK)]
ACR15WK_M <- pheno_M$ACR15WK[complete.cases(pheno_M$ACR15WK)]
ttest_ACR15WK <- t.test(ACR15WK_F, ACR15WK_M)
# 	Welch Two Sample t-test
#
# data:  pheno_F$ACR15WK and pheno_M$ACR15WK
# t = -8.469, df = 124.03, p-value = 5.972e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -3458.218 -2148.003
# sample estimates:
# mean of x mean of y
# 2356.728  5159.839
#
# standard error of mean
sem_ACR15WK_F<-sd(ACR15WK_F)/sqrt(length(ACR15WK_F)) # 153.3739
sem_ACR15WK_M<-sd(ACR15WK_M)/sqrt(length(ACR15WK_M)) # 293.303

