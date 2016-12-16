## Comparing Type IV Collagen expressions
# In this initial analysis we will be comparing Type IV Collagen a1 ,a2, a3, and a4
# We will be creating a corelation map between the following and using 8 colours to represent each founders:

#Col4a1 vs Col4a2
#Col4a1 vs Col4a3
#Col4a1 vs Col4a4

#Col4a2 vs Col4a3
#Col4a2 vs Col4a4
#Col4a3 vs Col4a4

#ENSEMBL ID
#Col4a1: ENSMUSG00000031502
#Col4a2: ENSMUSG00000031503
#Col4a3: ENSMUSG00000079465
#Col4a4: ENSMUSG00000067158
#Col4a5: ENSMUSG00000031274

#Prior to, on cadillac created a list of civet run list 
#ls civet_run/ > civet_list.txt

library(DOQTL)
library(ggplot2)
library(grid)
setwd("/hpcdata/ytakemon/Col4a5xDO")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
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
#keep only interested columns
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Alb10WK_log","Alb15WK_log","ACR6WK_log", "ACR10WK_log", "ACR15WK_log")]

#RankZ
RNA_seqZ <- t(RNA_seq)
RNA_seqZ <- apply(RNA_seqZ, 2, rankZ)
RNA_seqZ <- t(RNA_seqZ)

#sample
civet_dir <- list.files("./civet_run/", full.names = T)
civet_sample_list <- read.delim("./Sample_list/RNAseq_1415_sample_list.txt", header = FALSE) #formatted as civet 1415-0000
genoprob_names <- rownames(best.genoprobs.192)
Complete_set_names <- as.data.frame(civet_dir)
Complete_set_names$civet_names <- civet_sample_list
Complete_set_names$genoprob_names <- genoprob_names
colnames(Complete_set_names) <- names(Complete_set_names)
temp <- read.delim(file = paste(Complete_set_names[1, "civet_dir"],"/", "gbrs.quantified.diploid.genes.tpm", sep =""), header = TRUE, sep ="\t")
rownames(temp) <- temp$locus
temp$locus <- NULL
col_name <- colnames(temp)
locus_names <- rownames(temp)

#	get gene allele for each sample
#Gene_allele <- array(0, c(192, 46517), dimnames = list( Complete_set_names$genoprob_names, locus_names))
#for (i in 1:192){
#	temp <- read.delim(file = paste(Complete_set_names[i, "civet_dir"],"/", "gbrs.quantified.diploid.genes.tpm", sep =""), header = TRUE, sep ="\t")
#	Gene_allele[i,] <- as.character(temp$notes)
#}
#save(Gene_allele, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")
load(file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Gene_allele.Rdata")

#check to see Col4 genes are present
geneID <- colnames(RNA_seqZ)
geneID[geneID == "ENSMUSG00000031502"] #col1
geneID[geneID == "ENSMUSG00000031503"] #col2
geneID[geneID == "ENSMUSG00000079465"] #col3
geneID[geneID == "ENSMUSG00000067158"] #col4 
geneID[geneID == "ENSMUSG00000031274"] #col5
#Lookin good

#Extract Col4 genes into 1 data frame
Col_names <- c("Col4a1_tpm", "Col4a1_allele", "Col4a2_tpm", "Col4a2_allele", "Col4a3_tpm", "Col4a3_allele", "Col4a4_tpm", "Col4a4_allele", "Col4a5_tpm", "Col4a5_allele")
Col4_RNA_tpm <- array(0, c(192, 10), dimnames = list(Complete_set_names$genoprob_names, Col_names))
for (i in 1:192){
	#Col4a1
	Col4_RNA_tpm[i, "Col4a1_tpm"] <- RNA_seqZ[i, "ENSMUSG00000031502"]
	Col4_RNA_tpm[i, "Col4a1_allele"] <- Gene_allele[i, "ENSMUSG00000031502"]
	#Col4a2
	Col4_RNA_tpm[i, "Col4a2_tpm"] <- RNA_seqZ[i, "ENSMUSG00000031503"]
	Col4_RNA_tpm[i, "Col4a2_allele"] <- Gene_allele[i, "ENSMUSG00000031503"]
	#Col4a3
	Col4_RNA_tpm[i, "Col4a3_tpm"] <- RNA_seqZ[i, "ENSMUSG00000079465"]
	Col4_RNA_tpm[i, "Col4a3_allele"] <- Gene_allele[i, "ENSMUSG00000079465"]
	#Col4a4
	Col4_RNA_tpm[i, "Col4a4_tpm"] <- RNA_seqZ[i, "ENSMUSG00000067158"]
	Col4_RNA_tpm[i, "Col4a4_allele"] <- Gene_allele[i, "ENSMUSG00000067158"]
	#Col4a5
	Col4_RNA_tpm[i, "Col4a5_tpm"] <- RNA_seqZ[i, "ENSMUSG00000031274"]
	Col4_RNA_tpm[i, "Col4a5_allele"] <- Gene_allele[i, "ENSMUSG00000031274"]
}
Col4_RNA_tpm <- as.data.frame(Col4_RNA_tpm)

#	some a1 & a2 AND a3 & a4 are have different alleles...
#	a1 != a2 (3 samples)
#                 Col4a1_tpm Col4a1_allele       Col4a2_tpm Col4a2_allele
#X1415.0092 1.81306086486753            AB 1.58945962725693            BH
#X1415.0293 1.71389972794688            AB 1.48453464320424            BH
#X1415.1019 2.10789353593388            BH 1.90886670005156            BF

# a3 != a4 (only 2 samples)
#                 Col4a3_tpm Col4a3_allele       Col4a4_tpm Col4a4_allele
#X1415.0084 1.88445210164078            AB  1.6677723055656            BD
#X1415.0229 1.89085303106325            AB 1.74047984079122            BD

#1-3 samples that are not BX allele, so throw them out
#Clean Col4a1 allele
temp  <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "AB",]
temp2 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BB",]
temp3 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BC",]
temp4 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BD",]
temp5 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BE",]
temp6 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BF",]
temp7 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BG",]
temp8 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a1_allele == "BH",]
Col4_RNA_tpm_col1 <- rbind(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8) #dim 184  10
temp  <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "AB",]
temp2 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BB",]
temp3 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BC",]
temp4 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BD",]
temp5 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BE",]
temp6 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BF",]
temp7 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BG",]
temp8 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a2_allele == "BH",]
Col4_RNA_tpm_col2 <- rbind(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8) #dim 184  10
temp  <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "AB",]
temp2 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BB",]
temp3 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BC",]
temp4 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BD",]
temp5 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BE",]
temp6 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BF",]
temp7 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BG",]
temp8 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a3_allele == "BH",]
Col4_RNA_tpm_col3 <- rbind(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8) #dim 192  10
temp  <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "AB",]
temp2 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BB",]
temp3 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BC",]
temp4 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BD",]
temp5 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BE",]
temp6 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BF",]
temp7 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BG",]
temp8 <-  Col4_RNA_tpm[Col4_RNA_tpm$Col4a4_allele == "BH",]
Col4_RNA_tpm_col4 <- rbind(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8) #dim 192  10

#convert TPM values into numeric
col_dflist <- list(Col4_RNA_tpm_col1, Col4_RNA_tpm_col2, Col4_RNA_tpm_col3, Col4_RNA_tpm_col4)
col_dfnames <- c("Col4_RNA_tpm_col1", "Col4_RNA_tpm_col2", "Col4_RNA_tpm_col3", "Col4_RNA_tpm_col4")
for (i in 1:4){
	col_dflist[[i]]$Col4a1_tpm <- as.numeric(as.character(col_dflist[[i]]$Col4a1_tpm))
	col_dflist[[i]]$Col4a2_tpm <- as.numeric(as.character(col_dflist[[i]]$Col4a2_tpm))
	col_dflist[[i]]$Col4a3_tpm <- as.numeric(as.character(col_dflist[[i]]$Col4a3_tpm))
	col_dflist[[i]]$Col4a4_tpm <- as.numeric(as.character(col_dflist[[i]]$Col4a4_tpm))
	col_dflist[[i]]$Col4a5_tpm <- as.numeric(as.character(col_dflist[[i]]$Col4a5_tpm))
	assign(col_dfnames[i], col_dflist[[i]])
}

#Consolidate in to 1 data frame with dim 184 10
#order
Col4_RNA_tpm_col1 <- Col4_RNA_tpm_col1[ order(rownames(Col4_RNA_tpm_col1)),]
Col4_RNA_tpm_col2 <- Col4_RNA_tpm_col2[ order(rownames(Col4_RNA_tpm_col2)),]
Col4_RNA_tpm_col3 <- Col4_RNA_tpm_col3[ order(rownames(Col4_RNA_tpm_col3)),]
Col4_RNA_tpm_col4 <- Col4_RNA_tpm_col4[ order(rownames(Col4_RNA_tpm_col4)),]
#subset
Col4_RNA_tpm_col3 <- Col4_RNA_tpm_col3[rownames(Col4_RNA_tpm_col1),]
Col4_RNA_tpm_col4 <- Col4_RNA_tpm_col4[rownames(Col4_RNA_tpm_col1),]

#use Col4_RNA_tpm_col1 as reference
Col4a5_tpm_counts <- Col4_RNA_tpm_col1
save(Col4a5_tpm_counts, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Col4a5_rankZ_tpm_counts.Rdata")

#####################################################################
#Plot correlation between Type IV Collagen TPM rank counts:
#Col4a1 vs Col4a2
#Col4a1 vs Col4a3
#Col4a1 vs Col4a4

#Col4a2 vs Col4a3
#Col4a2 vs Col4a4
#Col4a3 vs Col4a4
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Col4a5_rankZ_tpm_counts.Rdata")
pheno <- pheno[rownames(Col4a5_tpm_counts),]
Col4a5_tpm_counts$Sex <- pheno$Sex
Col4a5_tpm_counts_F <- Col4a5_tpm_counts[Col4a5_tpm_counts$Sex == "F",]
Col4a5_tpm_counts_M <- Col4a5_tpm_counts[Col4a5_tpm_counts$Sex == "M",]

#Col4a1 vs Col4a2
fit <- lm(Col4a2_tpm ~ Col4a1_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ Col4a1_tpm, Col4a5_tpm_counts_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ Col4a1_tpm, Col4a5_tpm_counts_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot12 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a2_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , x = 1.9, y = 2.3, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a1 vs Col4a2 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot12_sex <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a2_tpm)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , x = 1.9, y = 2.3, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , x = 1.9, y = 2.25, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( colour = "Sex", title = "Col4a1 vs Col4a2 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

#png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a1.v.Col4a2.png", width = 1500, height = 1000, res = 100)
#print(ggplot12)
#dev.off()

#Col4a1 vs Col4a3
fit <- lm(Col4a3_tpm ~ Col4a1_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Col4a1_tpm, Col4a5_tpm_counts_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Col4a1_tpm, Col4a5_tpm_counts_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot13 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a3_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , x = 1.9, y = 2.1, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a1 vs Col4a3 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot13_alt <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a3_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , x = 1.9, y = 2.1, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a1 vs Col4a3 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot13_sex <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a3_tpm)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , x = 1.9, y = 2.1, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , x = 1.9, y = 2.05, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( colour = "Sex", title = "Col4a1 vs Col4a3 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

#png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a1.v.Col4a3.png", width = 1500, height = 1000, res = 100)
#print(ggplot13)
#dev.off()

#Col4a1 vs Col4a4
fit <- lm(Col4a4_tpm ~ Col4a1_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Col4a1_tpm, Col4a5_tpm_counts_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Col4a1_tpm, Col4a5_tpm_counts_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot14 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a4_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , x = 1.9, y = 1.9, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a1 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot14_alt <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a4_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , x = 1.9, y = 1.9, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a1 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot14_sex <- ggplot(Col4a5_tpm_counts, aes(x = Col4a1_tpm, y = Col4a4_tpm)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , x = 1.9, y = 1.9, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , x = 1.9, y = 1.85, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( colour = "Sex", title = "Col4a1 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 


#png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a1.v.Col4a4.png", width = 1500, height = 1000, res = 100)
#print(ggplot14)
#dev.off()

#Col4a2 vs Col4a3
fit <- lm(Col4a3_tpm ~ Col4a2_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Col4a2_tpm, Col4a5_tpm_counts_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Col4a2_tpm, Col4a5_tpm_counts_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot23 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a2_tpm, y = Col4a3_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , x = 1.7, y = 2.1, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a2 vs Col4a3 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot23_alt <- ggplot(Col4a5_tpm_counts, aes(x = Col4a2_tpm, y = Col4a3_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , x = 1.7, y = 2.1, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a2 vs Col4a3 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot23_sex <- ggplot(Col4a5_tpm_counts, aes(x = Col4a2_tpm, y = Col4a3_tpm)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , x = 1.7, y = 2.1, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , x = 1.7, y = 2.05, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( colour = "Sex", title = "Col4a2 vs Col4a3 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

#png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a2.v.Col4a3.png", width = 1500, height = 1000, res = 100)
#print(ggplot23)
#dev.off()

#Col4a2 vs Col4a4
fit <- lm(Col4a4_tpm ~ Col4a2_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Col4a2_tpm, Col4a5_tpm_counts_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Col4a2_tpm, Col4a5_tpm_counts_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot24 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a2_tpm, y = Col4a4_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , x = 1.7, y = 1.9, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a2 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot24_alt <- ggplot(Col4a5_tpm_counts, aes(x = Col4a2_tpm, y = Col4a4_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , x = 1.7, y = 1.9, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a2 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

ggplot24_sex <- ggplot(Col4a5_tpm_counts, aes(x = Col4a2_tpm, y = Col4a4_tpm)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , x = 1.7, y = 1.9, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , x = 1.7, y = 1.85, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( colour = "Sex", title = "Col4a2 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

#png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a2.v.Col4a4.png", width = 1500, height = 1000, res = 100)
#print(ggplot24)
#dev.off()

#Col4a3 vs Col4a4
fit <- lm(Col4a4_tpm ~ Col4a3_tpm, Col4a5_tpm_counts)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Col4a3_tpm, Col4a5_tpm_counts_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Col4a3_tpm, Col4a5_tpm_counts_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ggplot34 <- ggplot(Col4a5_tpm_counts, aes(x = Col4a3_tpm, y = Col4a4_tpm)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , x = 1.8, y = 1.9, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( title = "Col4a3 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ggplot34_sex <- ggplot(Col4a5_tpm_counts, aes(x = Col4a3_tpm, y = Col4a4_tpm)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , x = 1.8, y = 1.9, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , x = 1.8, y = 1.85, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_x_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	labs( colour = "Sex", title = "Col4a3 vs Col4a4 TPM") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a3.v.Col4a4.png", width = 1500, height = 1000, res = 100)
#print(ggplot34)
#dev.off()

library(grid)
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_allele1.png", width = 1500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(3, 3)))
print(ggplot12, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot13, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(ggplot14, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(ggplot23, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(ggplot24, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(ggplot34, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_allele2.png", width = 1500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(3, 3)))
print(ggplot12, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot13_alt, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(ggplot14_alt, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(ggplot23_alt, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(ggplot24_alt, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(ggplot34, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_sex.png", width = 1500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(3, 3)))
print(ggplot12_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot13_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(ggplot14_sex, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(ggplot23_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(ggplot24_sex, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(ggplot34_sex, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
dev.off()

#	Plot cor between Type IV collagen and phenotypes
png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Alb_hist.png", width = 1000, height = 1500, res = 100)
	layout(matrix(1:6,3,2))
	hist(pheno$Alb6WK)
	hist(pheno$Alb6WK_log)
	hist(pheno$Alb10WK)
	hist(pheno$Alb10WK_log)
	hist(pheno$Alb15WK)
	hist(pheno$Alb15WK_log)
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Creat_hist.png", width = 1000, height = 1500, res = 100)
	layout(matrix(1:6,3,2))
	hist(pheno$Creat6WK)
	hist(pheno$Creat6WK_log)
	hist(pheno$Creat10WK)
	hist(pheno$Creat10WK_log)
	hist(pheno$Creat15WK)
	hist(pheno$Creat15WK_log)
dev.off()

#subset pheno
pheno_sub <- pheno[rownames(Col4a5_tpm_counts),]
Col4a5_tmp_pheno <- Col4a5_tpm_counts
Col4a5_tmp_pheno$Sex <- pheno_sub$Sex
Col4a5_tmp_pheno$C2_log <- pheno_sub$C2_log
Col4a5_tmp_pheno$Alb6WK_log <- pheno_sub$Alb6WK_log
Col4a5_tmp_pheno$Alb10WK_log <- pheno_sub$Alb10WK_log
Col4a5_tmp_pheno$Alb15WK_log <- pheno_sub$Alb15WK_log
Col4a5_tmp_pheno$ACR6WK_log <- pheno_sub$ACR6WK_log
Col4a5_tmp_pheno$ACR10WK_log <- pheno_sub$ACR10WK_log
Col4a5_tmp_pheno$ACR15WK_log <- pheno_sub$ACR15WK_log
Col4a5_tmp_pheno$Alb6WK <- pheno_sub$Alb6WK
Col4a5_tmp_pheno$Alb10WK <- pheno_sub$Alb10WK
Col4a5_tmp_pheno$Alb15WK <- pheno_sub$Alb15WK
Col4a5_tmp_pheno$Creat6WK_log <- pheno_sub$Creat6WK_log
Col4a5_tmp_pheno$Creat10WK_log <- pheno_sub$Creat10WK_log
Col4a5_tmp_pheno$Creat15WK_log <- pheno_sub$Creat15WK_log
Col4a5_tmp_pheno$ACR6WK <- pheno_sub$ACR6WK
Col4a5_tmp_pheno$ACR10WK <- pheno_sub$ACR10WK
Col4a5_tmp_pheno$ACR15WK <- pheno_sub$ACR15WK

save(Col4a5_tmp_pheno, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Col4a5_rankZ_tmp_pheno.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/Col4a5_rankZ_tmp_pheno.Rdata")

#	Col4 TPM vs GFR
data <- Col4a5_tmp_pheno[ complete.cases( Col4a5_tmp_pheno$C2_log) ,]
data_F <- data[data$Sex == "F",]
data_M <- data[data$Sex == "M",]

#Col4a1 vs GFR
fit <- lm(Col4a1_tpm ~ C2_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ C2_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ C2_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")


gfr1 <- ggplot(data, aes(y = Col4a1_tpm, x = C2_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.6, x = 4.6, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( title = "Col4a1 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

gfr1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = C2_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.6, x = 4.6, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.55, x = 4.6, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a2 vs GFR
fit <- lm(Col4a2_tpm ~ C2_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ C2_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ C2_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

gfr2 <- ggplot(data, aes(y = Col4a2_tpm, x = C2_log)) +
	geom_smooth( method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.4, x = 4.6, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( title = "Col4a2 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

gfr2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = C2_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.4, x = 4.6, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.35, x = 4.6, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( colour = "Sex", title = "Col4a2 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a3 vs GFR
fit <- lm(Col4a3_tpm ~ C2_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ C2_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ C2_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

gfr3 <- ggplot(data, aes(y = Col4a3_tpm, x = C2_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.0, x = 4.6, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( title = "Col4a3 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

gfr3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = C2_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.0, x = 4.6, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.98, x = 4.6, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( colour = "Sex", title = "Col4a3 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a4 vs GFR
fit <- lm(Col4a4_tpm ~ C2_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ C2_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ C2_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

gfr4 <- ggplot(data, aes(y = Col4a4_tpm, x = C2_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.9, x = 4.6, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( title = "Col4a4 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

gfr4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = C2_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.9, x = 4.6, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.88, x = 4.6, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "log-transformed C2 GFR", breaks = seq(0, 7.5, by = 0.2)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs log C2 GFR") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_GFR_by_Allele.png", width = 2000, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 2)))
print(gfr1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(gfr2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gfr3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(gfr4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_GFR_by_sex.png", width = 2000, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 2)))
print(gfr1_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(gfr2_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gfr3_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(gfr4_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()


#	Col4 TPM vs Alb6WK
data <- Col4a5_tmp_pheno[ complete.cases( Col4a5_tmp_pheno$ACR6WK_log), ]
data_F <- data[data$Sex == "F",]
data_M <- data[data$Sex == "M",]


#Col4a1 vs ACR6
fit <- lm(Col4a1_tpm ~ ACR6WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a1_tpm ~ ACR6WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ ACR6WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR6_1 <- ggplot(data, aes(y = Col4a1_tpm, x = ACR6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.55, x = 4.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 4.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( title = "Col4a1 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR6_1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = ACR6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.6, x = 4.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.55, x = 4.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 4.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a1 vs Alb6
fit <- lm(Col4a1_tpm ~ Alb6WK_log + Creat6WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a1_tpm ~ Alb6WK_log + Creat6WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ Alb6WK_log + Creat6WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb6_1 <- ggplot(data, aes(y = Col4a1_tpm, x = Alb6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.55, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 1.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( title = "Col4a1 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb6_1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = Alb6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.6, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.55, x = 1.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 1.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a2 vs ACR6
fit <- lm(Col4a2_tpm ~ ACR6WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a2_tpm ~ ACR6WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ ACR6WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR6_2 <- ggplot(data, aes(y = Col4a2_tpm, x = ACR6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.35, x = 4.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 4.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( title = "Col4a2 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
ACR6_2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = ACR6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.4, x = 4.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.35, x = 4.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 4.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a2 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a2 vs Alb6
fit <- lm(Col4a2_tpm ~ Alb6WK_log + Creat6WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a2_tpm ~ Alb6WK_log + Creat6WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ Alb6WK_log + Creat6WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb6_2 <- ggplot(data, aes(y = Col4a2_tpm, x = Alb6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.35, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 1.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( title = "Col4a2 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb6_2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = Alb6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.4, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.35, x = 1.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 1.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a2 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a3 vs ACR6
fit <- lm(Col4a3_tpm ~ ACR6WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a3_tpm ~ ACR6WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ ACR6WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR6_3 <- ggplot(data, aes(y = Col4a3_tpm, x = ACR6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.12, x = 4.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 4.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( title = "Col4a3 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR6_3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = ACR6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.12, x = 4.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.10, x = 4.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.08, x = 4.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a3 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
#Col4a3 vs Alb6
fit <- lm(Col4a3_tpm ~ Alb6WK_log + Creat6WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a3_tpm ~ Alb6WK_log + Creat6WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Alb6WK_log + Creat6WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb6_3 <- ggplot(data, aes(y = Col4a3_tpm, x = Alb6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.12, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 1.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( title = "Col4a3 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb6_3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = Alb6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.12, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.10, x = 1.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.08, x = 1.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a3 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a4 vs ACR6
fit <- lm(Col4a4_tpm ~ ACR6WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a4_tpm ~ ACR6WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ ACR6WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR6_4 <- ggplot(data, aes(y = Col4a4_tpm, x = ACR6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.82, x = 4.8, label = fit_form, fontface = "bold", size = 3) +
	annotate( "text" , y = 1.8, x = 4.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( title = "Col4a4 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR6_4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = ACR6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.84, x = 4.8, label = fit_form, fontface = "bold", size = 3) +
	annotate( "text" , y = 1.82, x = 4.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.80, x = 4.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR6WK_log", breaks = seq(0, 8.8, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs ACR6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
#Col4a4 vs Alb6
fit <- lm(Col4a4_tpm ~ Alb6WK_log + Creat6WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a4_tpm ~ Alb6WK_log + Creat6WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Alb6WK_log + Creat6WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb6_4 <- ggplot(data, aes(y = Col4a4_tpm, x = Alb6WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.82, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.8, x = 1.8, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( title = "Col4a4 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb6_4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = Alb6WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.84, x = 1.8, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.82, x = 1.8, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.80, x = 1.8, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb6WK_log", breaks = seq(-1.0, 5.2, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs Alb6WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_Albumin6WK_log_by_allele.png", width = 2500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 4)))
print(ACR6_1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR6_2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ACR6_3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ACR6_4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(Alb6_1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Alb6_2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Alb6_3, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(Alb6_4, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_Albumin6WK_log_by_sex.png", width = 2500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 4)))
print(ACR6_1_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR6_2_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ACR6_3_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ACR6_4_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(Alb6_1_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Alb6_2_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Alb6_3_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(Alb6_4_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
dev.off()

#	Col4 TPM vs Alb10WK
data <- Col4a5_tmp_pheno[ complete.cases( Col4a5_tmp_pheno$ACR10WK_log), ]
data_F <- data[data$Sex == "F",]
data_M <- data[data$Sex == "M",]

#Col4a1 vs ACR6
fit <- lm(Col4a1_tpm ~ ACR10WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a1_tpm ~ ACR10WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ ACR10WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR10_1 <- ggplot(data, aes(y = Col4a1_tpm, x = ACR10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.55, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a1 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR10_1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = ACR10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.55, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.45, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a1 vs Alb10
fit <- lm(Col4a1_tpm ~ Alb10WK_log + Creat10WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a1_tpm ~ Alb10WK_log + Creat10WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ Alb10WK_log + Creat10WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb10_1 <- ggplot(data, aes(y = Col4a1_tpm, x = Alb10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.55, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a1 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb10_1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = Alb10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.55, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.45, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a2 vs ACR10
fit <- lm(Col4a2_tpm ~ ACR10WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a2_tpm ~ ACR10WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ ACR10WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR10_2 <- ggplot(data, aes(y = Col4a2_tpm, x = ACR10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.35, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a2 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR10_2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = ACR10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.35, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.25, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a2 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
#Col4a2 vs Alb10
fit <- lm(Col4a2_tpm ~ Alb10WK_log + Creat10WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a2_tpm ~ Alb10WK_log + Creat10WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ Alb10WK_log + Creat10WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb10_2 <- ggplot(data, aes(y = Col4a2_tpm, x = Alb10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.35, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a2 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb10_2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = Alb10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.35, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.25, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a2 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a3 vs ACR10
fit <- lm(Col4a3_tpm ~ ACR10WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a3_tpm ~ ACR10WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ ACR10WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR10_3 <- ggplot(data, aes(y = Col4a3_tpm, x = ACR10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.12, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a3 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
ACR10_3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = ACR10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.12, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.08, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a3 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	

#Col4a3 vs Alb10
fit <- lm(Col4a3_tpm ~ Alb10WK_log + Creat10WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a3_tpm ~ Alb10WK_log + Creat10WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Alb10WK_log + Creat10WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb10_3 <- ggplot(data, aes(y = Col4a3_tpm, x = Alb10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.12, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a3 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb10_3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = Alb10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.12, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.08, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a3 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a4 vs ACR10
fit <- lm(Col4a4_tpm ~ ACR10WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a4_tpm ~ ACR10WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ ACR10WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR10_4 <- ggplot(data, aes(y = Col4a4_tpm, x = ACR10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.82, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.8, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a4 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR10_4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = ACR10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.82, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.8, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.78, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR10WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs ACR10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
#Col4a4 vs Alb10
fit <- lm(Col4a4_tpm ~ Alb10WK_log + Creat10WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a4_tpm ~ Alb10WK_log + Creat10WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Alb10WK_log + Creat10WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb10_4 <- ggplot(data, aes(y = Col4a4_tpm, x = Alb10WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.82, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.8, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a4 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb10_4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = Alb10WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.82, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.8, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.78, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb10WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs Alb10WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_Albumin10WK_log_by_allele.png", width = 2500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 4)))
print(ACR10_1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR10_2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ACR10_3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ACR10_4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(Alb10_1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Alb10_2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Alb10_3, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(Alb10_4, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_Albumin10WK_log_by_sex.png", width = 2500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 4)))
print(ACR10_1_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR10_2_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ACR10_3_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ACR10_4_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(Alb10_1_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Alb10_2_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Alb10_3_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(Alb10_4_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
dev.off()

#	Col4 TPM vs Alb15WK
data <- Col4a5_tmp_pheno[ complete.cases( Col4a5_tmp_pheno$ACR15WK_log), ]
data_F <- data[data$Sex == "F",]
data_M <- data[data$Sex == "M",]

#Col4a1 vs ACR6
fit <- lm(Col4a1_tpm ~ ACR15WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a1_tpm ~ ACR15WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ ACR15WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR15_1 <- ggplot(data, aes(y = Col4a1_tpm, x = ACR15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.55, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a1 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

ACR15_1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = ACR15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.55, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.45, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a1 vs Alb15
fit <- lm(Col4a1_tpm ~ Alb15WK_log + Creat15WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a1_tpm ~ Alb15WK_log + Creat15WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a1_tpm ~ Alb15WK_log + Creat15WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb15_1 <- ggplot(data, aes(y = Col4a1_tpm, x = Alb15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a1_allele), size = 2) +
	annotate( "text" , y = 2.55, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a1 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb15_1_sex <- ggplot(data, aes(y = Col4a1_tpm, x = Alb15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.55, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.5, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.45, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a1 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a1 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a2 vs ACR15
fit <- lm(Col4a2_tpm ~ ACR15WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a2_tpm ~ ACR15WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ ACR15WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR15_2 <- ggplot(data, aes(y = Col4a2_tpm, x = ACR15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.35, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a2 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
ACR15_2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = ACR15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.35, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.25, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a2 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a2 vs Alb15
fit <- lm(Col4a2_tpm ~ Alb15WK_log + Creat15WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a2_tpm ~ Alb15WK_log + Creat15WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a2_tpm ~ Alb15WK_log + Creat15WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb15_2 <- ggplot(data, aes(y = Col4a2_tpm, x = Alb15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a2_allele), size = 2) +
	annotate( "text" , y = 2.35, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a2 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb15_2_sex <- ggplot(data, aes(y = Col4a2_tpm, x = Alb15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.35, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.3, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.25, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a2 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", itle = "Col4a2 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  


#Col4a3 vs ACR15
fit <- lm(Col4a3_tpm ~ ACR15WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a3_tpm ~ ACR15WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ ACR15WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR15_3 <- ggplot(data, aes(y = Col4a3_tpm, x = ACR15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.12, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a3 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
ACR15_3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = ACR15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.12, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.08, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a3 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a3 vs Alb15
fit <- lm(Col4a3_tpm ~ Alb15WK_log + Creat15WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a3_tpm ~ Alb15WK_log + Creat15WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a3_tpm ~ Alb15WK_log + Creat15WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb15_3 <- ggplot(data, aes(y = Col4a3_tpm, x = Alb15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a3_allele), size = 2) +
	annotate( "text" , y = 2.12, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a3 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb15_3_sex <- ggplot(data, aes(y = Col4a3_tpm, x = Alb15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 2.12, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.1, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 2.08, x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a3 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a3 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a4 vs ACR15
fit <- lm(Col4a4_tpm ~ ACR15WK_log , data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a4_tpm ~ ACR15WK_log , data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ ACR15WK_log , data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

ACR15_4 <- ggplot(data, aes(y = Col4a4_tpm, x = ACR15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.85, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.83, x = 6.0, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( title = "Col4a4 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  
	
ACR15_4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = ACR15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.85, x = 6.0, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.83, x = 6.0, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.81, x = 6.0, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "ACR15WK_log", breaks = seq(0, 9.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs ACR15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

#Col4a4 vs Alb15
fit <- lm(Col4a4_tpm ~ Alb15WK_log + Creat15WK_log, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
fit_form <- paste(formula(fit)[2], " ", formula(fit)[1], " ", formula(fit)[3], sep = "")

fit <- lm(Col4a4_tpm ~ Alb15WK_log + Creat15WK_log, data_F)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_F <- paste("Female = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

fit <- lm(Col4a4_tpm ~ Alb15WK_log + Creat15WK_log, data_M)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq_M <- paste("Male = ", slope,"x ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

Alb15_4 <- ggplot(data, aes(y = Col4a4_tpm, x = Alb15WK_log)) +
	geom_smooth(method = lm) +
	geom_point( aes(colour = Col4a4_allele), size = 2) +
	annotate( "text" , y = 1.85, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.83, x = 2.5, label = eq, fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( title = "Col4a4 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

Alb15_4_sex <- ggplot(data, aes(y = Col4a4_tpm, x = Alb15WK_log)) +
	geom_smooth(aes(colour = factor(Sex)), method = lm) +
	geom_point( aes(colour = Sex), size = 2) +
	annotate( "text" , y = 1.85, x = 2.5, label = fit_form, fontface = "bold", size = 3) + 
	annotate( "text" , y = 1.83, x = 2.5, label = eq_F, colour = "#F8766D", fontface = "bold", size = 3) + 
	annotate( "text" , y = 1., x = 2.5, label = eq_M, colour = "#00BFC4", fontface = "bold", size = 3) + 
	guides( colour = "legend") +
	scale_y_continuous( "Col4a4 TPM", breaks = seq(0, 2.5, by = 0.1)) +
	scale_x_continuous( "Alb15WK_log", breaks = seq(-1.0, 5.5, by = 0.5)) +
	labs( colour = "Sex", title = "Col4a4 TPM vs Alb15WK_log") +
	theme(legend.position = "right", plot.title = element_text(hjust = 0.5))  

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_Albumin15WK_log_by_allele.png", width = 2500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 4)))
print(ACR15_1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR15_2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ACR15_3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ACR15_4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(Alb15_1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Alb15_2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Alb15_3, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(Alb15_4, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
dev.off()

png("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4aX_Albumin15WK_log_by_sex.png", width = 2500, height = 1000, res = 100)
pushViewport(viewport(layout = grid.layout(2, 4)))
print(ACR15_1_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ACR15_2_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ACR15_3_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ACR15_4_sex, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(Alb15_1_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Alb15_2_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Alb15_3_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(Alb15_4_sex, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
dev.off()
##################
