# Run this script
# R-3.6.0 --no-save < Phenotype_RNA_cor_bySex_M.R > Rout/Phenotype_RNA_cor_bySex_M.Rout 2>&1 &

# Create a heatmap of correlation between gene expression and ACR @ 15wks and GFR at 14wks

library(pacman)
p_load(tidyverse, broom)

cor_dir <- "/projects/marralab/ytakemon_prj/Col4a5/"
load(paste0(cor_dir,"Data/RNAseq/RNA_seq_rankZ_tpm.Rdata"))
load(paste0(cor_dir,"Data/RNAseq/RNA_seq_tpm.Rdata"))
load(paste0(cor_dir,"Data/Col4a5xDO_192data_YT.Rdata"))

# get a list of genes that have no expression and remove
zero <- NULL
for(i in 1:ncol(RNA_seq)){
  gene <- colnames(RNA_seq)[i]
  if(all(RNA_seq[,gene] == 0)){
    zero <- c(zero, gene)
  } else {
    next
  }
}

select_pheno <- Pheno %>% mutate(ACR15WK_log = log(ACR15WK)) %>%
  filter(Sex == "M") %>%
  select(MouseID, ACR15WK_log, C2_log)
select_RNA_seqZ <- RNA_seqZ %>% as_tibble(., rownames = "MouseID") %>%
    filter(MouseID %in% select_pheno$MouseID) %>%
    select(-all_of(zero)) #get rid of zeros

identical(select_RNA_seqZ$MouseID, select_pheno$MouseID) # good to compare directly

# Create function to get results
getRNACor <- function(gene, pheno, res){

  rna <- select_RNA_seqZ %>% select(MouseID, all_of(gene))
  df <- full_join(select_pheno, rna, by = "MouseID")

  fit <- cor.test(df[,pheno], df[,gene], use = "complete.obs", method = "pearson" ) %>%tidy()

  output <- fit %>% pull(res)
  return(output)
}

cor <- tibble(EnsID = colnames(select_RNA_seqZ)[-1],
              acr_cor = NA,
              acr_pval  = NA,
              gfr = NA,
              gfr_pval = NA) %>%
        mutate(
          acr_cor = map_dbl(EnsID, getRNACor, "ACR15WK_log", "estimate"),
          acr_pval = map_dbl(EnsID, getRNACor, "ACR15WK_log", "p.value"),
          gfr = map_dbl(EnsID, getRNACor, "C2_log", "estimate"),
          gfr_pval = map_dbl(EnsID, getRNACor, "C2_log", "p.value"))

write_csv(cor, paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk_M.csv"))

# ANALYSIS -------------------------------------------------------
cor <- read_csv(paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk_M.csv")) %>%
  mutate(acr_fdr = p.adjust(acr_pval, n = nrow(cor), method = "BH"),
         gfr_fdr = p.adjust(gfr_pval, n = nrow(cor), method = "BH")) %>%
  rename(gfr_cor = gfr)

write_csv(cor, paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk_M_wFDR.csv"))

cor_sig_both <- cor %>% filter(acr_fdr < 0.05 & gfr_fdr < 0.05)

cor %>%
  rename(ACR = acr_cor,
         GFR = gfr_cor) %>%
  select(EnsID, ACR, GFR) %>%
  pivot_longer(-EnsID, names_to = "Phenotype", values_to = "Cor") %>%
  ggplot(., aes(x = Cor, colour = Phenotype))+
    geom_freqpoly(bins = 20)+
    labs(x = "Correlation",
         y = "Frequency")

cor %>%
  rename(ACR = acr_cor,
         GFR = gfr_cor) %>%
  filter(acr_fdr < 0.05 & gfr_fdr < 0.05) %>%
  select(EnsID, ACR, GFR) %>%
  pivot_longer(-EnsID, names_to = "Phenotype", values_to = "Cor") %>%
  ggplot(., aes(x = Cor, colour = Phenotype))+
    geom_freqpoly(bins = 20)+
    labs(x = "Correlation",
         y = "Frequency")

pdf(paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_scatter_M.pdf"), height = 5, width = 7)
cor %>%
  rename(ACR = acr_cor,
         GFR = gfr_cor) %>%
  mutate(both_sig = acr_fdr < 0.05 & gfr_fdr < 0.05) %>%
  ggplot(., aes(x = ACR, y = GFR, colour = both_sig))+
    geom_point(alpha = 0.2) +
    scale_colour_manual(values = c("grey","dark red")) +
    labs(colour = "Significant")
dev.off()

# Create heatmap
p_load(pheatmap)
heat_df <- cor %>%
  rename(ACR = acr_cor,
         GFR = gfr_cor) %>%
  filter(acr_fdr < 0.05 & gfr_fdr < 0.05) %>%
  select(EnsID, ACR, GFR)
write_csv(heat_df, paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_heatmap_sig_M.csv"))

heat_mat_M <- heat_df %>% select(ACR, GFR) %>% as.matrix()
rownames(heat_mat_M) <- heat_df %>% pull(EnsID)

pdf(paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_heatmap_sig_M.pdf"), height = 5, width = 3)
pheatmap(heat_mat_M, cutree_rows = 2, scale = "none",
  annotation_names_row = FALSE,
  show_rownames = FALSE)
dev.off()

# Show F and M side by side
cor_M <- read_csv(paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk_M_wFDR.csv")) %>%
  filter(acr_fdr < 0.05 & gfr_fdr < 0.05) %>%
  select(EnsID, acr_cor, gfr_cor) %>%
  rename(acr_cor_M = acr_cor,
         gfr_cor_M = gfr_cor)

cor_F <- read_csv(paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk_F_wFDR.csv")) %>%
  filter(acr_fdr < 0.05 & gfr_fdr < 0.05) %>%
  select(EnsID, acr_cor, gfr_cor) %>%
  rename(acr_cor_F = acr_cor,
         gfr_cor_F = gfr_cor)

cor_both <- full_join(cor_M, cor_F) %>% filter(complete.cases(.))
write_csv(cor_both, paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_heatmap_sig_MF.csv"))
cor_both_mat <- cor_both %>% select(-EnsID) %>% as.matrix
rownames(cor_both_mat) <- cor_both$EnsID

pdf(paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_heatmap_sig_MF.pdf"), height = 5, width = 3)
pheatmap(cor_both_mat, cutree_rows = 2, scale = "none",
  method = "complete",
  annotation_names_row = FALSE,
  show_rownames = FALSE,
  na_col = "grey")
dev.off()
