# Yuka Takemon
# 08/29/18
# Consolidate main data in to one Rdata
setwd("/projects/korstanje-lab/ytakemon/Col4a5xDO/")

# Load individual data for consolidation
load("best.genoprobs.192.Rdata")
load("GM_snps.Rdata")
load("k.best.probs192.Rdata")
load("sex.covar.Rdata")
load("./RNAseq/RNA_seq_rankZ_tpm.Rdata")
load("./RNAseq/RNA_seq_tpm.Rdata")

# clean and gather important phenotype information
master_pheno <- read.delim("./Phenotype/1415_master_pheno.txt",
  sep = "\t",
  stringsAsFactors = FALSE)
small_pheno <- read.delim("./Phenotype/Minimal_shiny_pheno.txt",
  sep = "\t",
  stringsAsFactors = FALSE)

Pheno <- semi_join(master_pheno, small_pheno) %>%
  select(MouseID, Sex, DOB, SireGeneration, Weight,
    Alb6WK, Alb10WK, Alb15WK,
    Creat6WK, Creat10WK, Creat15WK,
    C2, GFRDate)

Missed <- small_pheno[!(small_pheno$MouseID %in% Pheno$MouseID),]$MouseID
Missed <- master_pheno[master_pheno$MouseID %in% Missed,] %>%
  select(MouseID, Sex, DOB, SireGeneration, Weight,
    Alb6WK, Alb10WK, Alb15WK,
    Creat6WK, Creat10WK, Creat15WK,
    C2, GFRDate)

Pheno <- Pheno %>% bind_rows(Missed) %>%
  arrange(MouseID) %>%
  mutate(MouseID = rownames(best.genoprobs.192),
    ACR6WK = Alb6WK/Creat6WK * 1000,
    ACR10WK = Alb10WK/Creat10WK * 1000,
    ACR15WK = Alb15WK/Creat15WK * 1000,
    C2_log = log(C2))
rownames(Pheno) <- Pheno$Mouse.ID

# Add covar
Covar <- sex.covar %>%
  mutate(Generation = Pheno$SireGeneration)
rownames(Covar) <- Pheno$Mouse.ID)

# RNAseq data
mRNAexpr <- as.data.frame(RNA_seqZ)
mRNAraw <- as.data.frame(RNA_seq)


# Rename data for consolidation
Genoprobs <- best.genoprobs.192
Snps <- GM_snps
K <- K.probs
Pheno <- Pheno
Covar <- Covar

# Save as one big Rdata
save(Genoprobs, Snps, K, Pheno, Covar, mRNAraw, mRNAexpr, file = "Col4a5xDO_192data_YT.Rdata")

# Save all QTL in one ----------------------------------------------------------
setwd("~/Dropbox/Col4a5/Data/")
load("./QTL/qtl.GFR.log.C2.192.Rdata")
load("./QTL/qtl.log.Alb6WK.192.Rdata")
load("./QTL/qtl.log.Alb10WK.192.Rdata")
load("./QTL/qtl.log.Alb15WK.192.Rdata")
save(qtl.GFR.log.C2.192, qtl.log.Alb6WK.192, qtl.log.Alb10WK.192, qtl.log.Alb15WK.192,
    file = "Col4a5xDO_QTL_YT.Rdata")

# Save all QTL perms in one ----------------------------------------------------
load("./QTL/perms.1000.qtl.GFR.log.C2.192.Rdata")
load("./QTL/perms.1000.qtl.log.Alb6WK.192.Rdata")
load("./QTL/perms.1000.qtl.log.Alb10WK.192.Rdata")
load("./QTL/perms.1000.qtl.log.Alb15WK.192.Rdata")
save(perms.1000.qtl.GFR.log.C2.192, perms.1000.qtl.log.Alb10WK.192, perms.1000.qtl.log.Alb15WK.192, perms.1000.qtl.log.Alb6WK.192,
  file = "Col4a5xDO_QTLperm1000_YT.Rdata")

# Save all Gwas in one ---------------------------------------------------------
load("./GWAS/gwas.log.C2.Rdata")
load("./GWAS/gwas.Alb6WK.Rdata")
load("./GWAS/gwas.Alb10WK.Rdata")
load("./GWAS/gwas.Alb15WK.Rdata")
save(GWAS.Alb6WK, GWAS.Alb10WK, GWAS.Alb15WK, GWAS.log.C2,
  file = "Col4a5xDO_GWAS_YT.Rdata")

# Save all Gwas perm in one ----------------------------------------------------
GWAS.Alb6WK.perm <- read.delim("./GWAS/all.perm.Alb6.txt", sep = "\t", header = F)
GWAS.Alb10WK.perm <- read.delim("./GWAS/all.perm.Alb10.txt", sep = "\t")
GWAS.Alb15WK.perm <- read.delim("./GWAS/all.perm.Alb15.txt", sep = "\t")
GWAS.GFR.perm <- read.delim("./GWAS/all.perm.GFR.txt", sep = "\t")
names(GWAS.Alb6WK.perm) <- names(GWAS.Alb10WK.perm) <- names(GWAS.Alb15WK.perm) <- names(GWAS.GFR.perm) <- "perm1000"
save(GWAS.Alb6WK.perm, GWAS.Alb10WK.perm, GWAS.Alb15WK.perm, GWAS.GFR.perm,
  file = "Col4a5xDO_GWASperm1000_YT.Rdata")
