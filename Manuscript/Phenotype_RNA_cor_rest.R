# Run this script
# R-3.6.0 --no-save < Phenotype_RNA_cor_rest.R > Rout/Phenotype_RNA_cor_rest.Rout 2>&1 &

# Create a heatmap of correlation between gene expression and ACR @ 15wks and GFR at 14wks

library(pacman)
p_load(tidyverse, broom)

cor_dir <- "/projects/marralab/ytakemon_prj/Col4a5/"
load(paste0(cor_dir,"Data/RNAseq/RNA_seq_rankZ_tpm.Rdata"))
load(paste0(cor_dir,"Data/Col4a5xDO_192data_YT.Rdata"))

RNA_seqZ <- RNA_seqZ %>% as_tibble(., rownames = "MouseID")
select_pheno <- Pheno %>%
  mutate(ACR6WK_log = log(ACR6WK),
         ACR10WK_log = log(ACR10WK),
         ACR15WK_log = log(ACR15WK)) %>%
  select(MouseID, ACR6WK_log, ACR10WK_log, ACR15WK_log, C2_log)

identical(RNA_seqZ$MouseID, select_pheno$MouseID) # good to compare directly

# Create function to get results
getRNACor <- function(gene, pheno, res){

  rna <- RNA_seqZ %>% select(MouseID, all_of(gene))
  df <- full_join(select_pheno, rna, by = "MouseID")

  fit <- cor.test(df[,pheno], df[,gene], use = "complete.obs", method = "pearson" ) %>%tidy()

  output <- fit %>% pull(res)
  return(output)
}

cor <- tibble(EnsID = colnames(RNA_seqZ)[-1],
              acr6_cor = NA,
              acr6_pval  = NA,
              acr10_cor = NA,
              acr10_pval  = NA) %>%
        mutate(
          acr6_cor = map_dbl(EnsID, getRNACor, "ACR6WK_log", "estimate"),
          acr6_pval = map_dbl(EnsID, getRNACor, "ACR6WK_log", "p.value"),
          acr10_cor = map_dbl(EnsID, getRNACor, "ACR10WK_log", "estimate"),
          acr10_pval = map_dbl(EnsID, getRNACor, "ACR10WK_log", "p.value"))

write_csv(cor, paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_other_wk.csv"))
