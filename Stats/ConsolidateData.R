setwd("/projects/korstanje-lab/ytakemon/Col4a5xDO/")

# Load individual data for consolidation
load("best.genoprobs.192.Rdata")
load("GM_snps.Rdata")
load("k.best.probs192.Rdata")
load("sex.covar.Rdata")

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
    ACR15WK = Alb15WK/Creat15WK * 1000)

# Add covar
Covar <- sex.covar %>%
  mutate(Generation = Pheno$SireGeneration)

# Rename data for consolidation
Genoprobs <- best.genoprobs.192
Snps <- GM_snps
K <- K.probs
Pheno <- Pheno
Covar <- Covar

# Save as one big Rdata
save(Genoprobs, Snps, K, Pheno, Covar, file = "Col4a5xDO_192data_YT.Rdata")
