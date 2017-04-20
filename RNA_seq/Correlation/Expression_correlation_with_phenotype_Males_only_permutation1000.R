setwd("/hpcdata/ytakemon/Col4a5xDO")

#Calculate correlation between phenotype and gene expression via tpm
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
RNA_seq <- as.data.frame(RNA_seq)
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
pheno <- pheno[,c("MouseID", "Sex", "C2_log", "Alb6WK_log","Creat6WK_log",
                  "Alb10WK_log","Creat10WK_log","Alb15WK_log","Creat15WK_log",
                  "ACR6WK_log", "ACR10WK_log", "ACR15WK_log")]
# Subset out males
pheno_M <- pheno[pheno$Sex == "M",] #95x12

#GFR
G_pheno <- pheno_M[complete.cases(pheno_M$C2_log),] #remove NAs
G_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(G_pheno),]
#randomize GFR
G_RNA_perm <- data.frame(perm_count = 1:1000, stringsAsFactors = FALSE)
G_RNA_perm$min_pval <- 0

for (i in 1:1000){
  #create temporary random data
  temp_G_pheno <- G_pheno
  temp_G_RNA_seq <- G_RNA_seq
  random <- sample(rownames(G_pheno)) #Take randome sampling of G_pheno's rownames
  rownames(temp_G_pheno) <- random
  temp_G_pheno <- temp_G_pheno[order(rownames(temp_G_pheno)),]
  #setup blank dataframes to hold data
  temp_RNA_GFR_cor <- array(0, c(length(colnames(temp_G_RNA_seq)),2),
                dimnames = list (colnames(temp_G_RNA_seq), c("GFR_correlation", "pval")))
  temp_RNA_GFR_cor <- as.data.frame(temp_RNA_GFR_cor)
  #calculate correlation for random set of data
  for (n in 1:length(colnames(temp_G_RNA_seq))){
    temp <- cor.test(temp_G_pheno$C2_log, temp_G_RNA_seq[,n])
    temp_RNA_GFR_cor[n, "GFR_correlation"] <- temp$estimate[[1]]
    temp_RNA_GFR_cor[n, "pval"] <- temp$p.value[[1]]
  }
  #find smallest pvalue for random set i
  G_RNA_perm[i, "min_pval"] <- min(temp_RNA_GFR_cor$pval, na.rm = TRUE)
  print(paste0("done with perm ", i, "."))
}
save(G_RNA_perm, file = "GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/G_RNA_perm.Rdata")

#ACR6
ACR6_pheno <- pheno_M[complete.cases(pheno_M$ACR6WK_log),] #remove NAs
ACR6_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(ACR6_pheno),]
#randomize GFR
ACR6_RNA_perm <- data.frame(perm_count = 1:1000, stringsAsFactors = FALSE)
ACR6_RNA_perm$min_pval <- 0

for (i in 1:1000){
  #create temporary random data
  temp_ACR6_pheno <- ACR6_pheno
  temp_ACR6_RNA_seq <- ACR6_RNA_seq
  random <- sample(rownames(ACR6_pheno)) #Take randome sampling of ACR6_pheno's rownames
  rownames(temp_ACR6_pheno) <- random
  temp_ACR6_pheno <- temp_ACR6_pheno[order(rownames(temp_ACR6_pheno)),]
  #setup blank dataframes to hold data
  temp_RNA_ACR6_cor <- array(0, c(length(colnames(temp_ACR6_RNA_seq)),2),
                dimnames = list (colnames(temp_ACR6_RNA_seq), c("ACR6_correlation", "pval")))
  temp_RNA_ACR6_cor <- as.data.frame(temp_RNA_ACR6_cor)
  #calculate correlation for random set of data
  for (n in 1:length(colnames(temp_ACR6_RNA_seq))){
    temp <- cor.test(temp_ACR6_pheno$ACR6WK_log, temp_ACR6_RNA_seq[,n])
    temp_RNA_ACR6_cor[n, "ACR6_correlation"] <- temp$estimate[[1]]
    temp_RNA_ACR6_cor[n, "pval"] <- temp$p.value[[1]]
  }
  #find smallest pvalue for random set i
  ACR6_RNA_perm[i, "min_pval"] <- min(temp_RNA_ACR6_cor$pval, na.rm = TRUE)
  print(paste0("done with perm ", i, "."))
}
save(ACR6_RNA_perm, file = "GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/ACR6_RNA_perm.Rdata")

#ACR10
ACR10_pheno <- pheno_M[complete.cases(pheno_M$ACR10WK_log),] #remove NAs
ACR10_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(ACR10_pheno),]
#randomize GFR
ACR10_RNA_perm <- data.frame(perm_count = 1:1000, stringsAsFactors = FALSE)
ACR10_RNA_perm$min_pval <- 0

for (i in 1:1000){
  #create temporary random data
  temp_ACR10_pheno <- ACR10_pheno
  temp_ACR10_RNA_seq <- ACR10_RNA_seq
  random <- sample(rownames(ACR10_pheno)) #Take randome sampling of ACR10_pheno's rownames
  rownames(temp_ACR10_pheno) <- random
  temp_ACR10_pheno <- temp_ACR10_pheno[order(rownames(temp_ACR10_pheno)),]
  #setup blank dataframes to hold data
  temp_RNA_ACR10_cor <- array(0, c(length(colnames(temp_ACR10_RNA_seq)),2),
                dimnames = list (colnames(temp_ACR10_RNA_seq), c("ACR10_correlation", "pval")))
  temp_RNA_ACR10_cor <- as.data.frame(temp_RNA_ACR10_cor)
  #calculate correlation for random set of data
  for (n in 1:length(colnames(temp_ACR10_RNA_seq))){
    temp <- cor.test(temp_ACR10_pheno$ACR10WK_log, temp_ACR10_RNA_seq[,n])
    temp_RNA_ACR10_cor[n, "ACR10_correlation"] <- temp$estimate[[1]]
    temp_RNA_ACR10_cor[n, "pval"] <- temp$p.value[[1]]
  }
  #find smallest pvalue for random set i
  ACR10_RNA_perm[i, "min_pval"] <- min(temp_RNA_ACR10_cor$pval, na.rm = TRUE)
  print(paste0("done with perm ", i, "."))
}
save(ACR10_RNA_perm, file = "GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/ACR10_RNA_perm.Rdata")

#ACR15
ACR15_pheno <- pheno_M[complete.cases(pheno_M$ACR15WK_log),] #remove NAs
ACR15_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(ACR15_pheno),]
#randomize GFR
ACR15_RNA_perm <- data.frame(perm_count = 1:1000, stringsAsFactors = FALSE)
ACR15_RNA_perm$min_pval <- 0

for (i in 1:1000){
  #create temporary random data
  temp_ACR15_pheno <- ACR15_pheno
  temp_ACR15_RNA_seq <- ACR15_RNA_seq
  random <- sample(rownames(ACR15_pheno)) #Take randome sampling of ACR15_pheno's rownames
  rownames(temp_ACR15_pheno) <- random
  temp_ACR15_pheno <- temp_ACR15_pheno[order(rownames(temp_ACR15_pheno)),]
  #setup blank dataframes to hold data
  temp_RNA_ACR15_cor <- array(0, c(length(colnames(temp_ACR15_RNA_seq)),2),
                dimnames = list (colnames(temp_ACR15_RNA_seq), c("ACR15_correlation", "pval")))
  temp_RNA_ACR15_cor <- as.data.frame(temp_RNA_ACR15_cor)
  #calculate correlation for random set of data
  for (n in 1:length(colnames(temp_ACR15_RNA_seq))){
    temp <- cor.test(temp_ACR15_pheno$ACR15WK_log, temp_ACR15_RNA_seq[,n])
    temp_RNA_ACR15_cor[n, "ACR15_correlation"] <- temp$estimate[[1]]
    temp_RNA_ACR15_cor[n, "pval"] <- temp$p.value[[1]]
  }
  #find smallest pvalue for random set i
  ACR15_RNA_perm[i, "min_pval"] <- min(temp_RNA_ACR15_cor$pval, na.rm = TRUE)
  print(paste0("done with perm ", i, "."))
}
save(ACR15_RNA_perm, file = "GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/ACR15_RNA_perm.Rdata")

#delta_ACR15_6
delta_ACR15_6_pheno <- pheno_M[complete.cases(pheno_M$delta_ACR15_6),]
delta_ACR15_6_RNA_seq <- RNA_seq[rownames(RNA_seq) %in% rownames(delta_ACR15_6_pheno),]
#randomize GFR
delta_ACR15_6_RNA_perm <- data.frame(perm_count = 1:1000, stringsAsFactors = FALSE)
delta_ACR15_6_RNA_perm$min_pval <- 0

for (i in 1:1000){
  #create temporary random data
  temp_delta_ACR15_6_pheno <- delta_ACR15_6_pheno
  temp_delta_ACR15_6_RNA_seq <- delta_ACR15_6_RNA_seq
  random <- sample(rownames(delta_ACR15_6_pheno)) #Take randome sampling of delta_ACR15_6_pheno's rownames
  rownames(temp_delta_ACR15_6_pheno) <- random
  temp_delta_ACR15_6_pheno <- temp_delta_ACR15_6_pheno[order(rownames(temp_delta_ACR15_6_pheno)),]
  #setup blank dataframes to hold data
  temp_RNA_delta_ACR15_6_cor <- array(0, c(length(colnames(temp_delta_ACR15_6_RNA_seq)),2),
                dimnames = list (colnames(temp_delta_ACR15_6_RNA_seq), c("delta_ACR15_6_correlation", "pval")))
  temp_RNA_delta_ACR15_6_cor <- as.data.frame(temp_RNA_delta_ACR15_6_cor)
  #calculate correlation for random set of data
  for (n in 1:length(colnames(temp_delta_ACR15_6_RNA_seq))){
    temp <- cor.test(temp_delta_ACR15_6_pheno$delta_ACR15_6, temp_delta_ACR15_6_RNA_seq[,n])
    temp_RNA_delta_ACR15_6_cor[n, "delta_ACR15_6_correlation"] <- temp$estimate[[1]]
    temp_RNA_delta_ACR15_6_cor[n, "pval"] <- temp$p.value[[1]]
    print(paste0("done with ", i, " of ", n, "."))
  }
  #find smallest pvalue for random set i
  delta_ACR15_6_RNA_perm[i, "min_pval"] <- min(temp_RNA_delta_ACR15_6_cor$pval, na.rm = TRUE)
  print(paste0("done with perm ", i, "."))
}
save(delta_ACR15_6_RNA_perm, file = "GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/delta_ACR15_6_RNA_perm.Rdata")
