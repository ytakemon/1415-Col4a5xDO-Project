# Run this script
# R-3.6.0 --no-save < Phenotype_RNA_cor.R > Rout/Phenotype_RNA_cor.Rout 2>&1 &

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
  print(i)
  gene <- colnames(RNA_seq)[i]
  if(all(RNA_seq[,gene] == 0)){
    zero <- c(zero, gene)
  } else {
    next
  }
}

RNA_seqZ <- RNA_seqZ %>% as_tibble(., rownames = "MouseID")
select_pheno <- Pheno %>% mutate(ACR15WK_log = log(ACR15WK)) %>%
  select(MouseID, ACR15WK_log, C2_log)

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
              acr_cor = NA,
              acr_pval  = NA,
              gfr = NA,
              gfr_pval = NA) %>%
        mutate(
          acr_cor = map_dbl(EnsID, getRNACor, "ACR15WK_log", "estimate"),
          acr_pval = map_dbl(EnsID, getRNACor, "ACR15WK_log", "p.value"),
          gfr = map_dbl(EnsID, getRNACor, "C2_log", "estimate"),
          gfr_pval = map_dbl(EnsID, getRNACor, "C2_log", "p.value"))

write_csv(cor, paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk.csv"))

cor <- read_csv(paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk.csv")) %>%
  mutate(acr_fdr = p.adjust(acr_pval, n = nrow(cor), method = "BH"),
         gfr_fdr = p.adjust(gfr_pval, n = nrow(cor), method = "BH")) %>%
  rename(gfr_cor = gfr_cor) %>%
  filter(!(EnsID %in% zero))

write_csv(cor, paste0(cor_dir,"Data/RNAseq/RNAseq_phenotype_cor_15wk_wFDR.csv"))

write_csv(cor, "~/Desktop/RNAseq_phenotype_cor_15wk_wFDR.csv")


cor_sig_both <- cor %>% filter(acr_fdr < 0.05 & gfr_fdr < 0.05)

pdf(paste0(cor_dir,"Results/RNAseq_cor/Gene_cor_all.pdf"), height = 3, width = 5)
cor %>%
  rename(ACR = acr_cor,
         GFR = gfr_cor) %>%
  select(EnsID, ACR, GFR) %>%
  pivot_longer(-EnsID, names_to = "Phenotype", values_to = "Cor") %>%
  ggplot(., aes(x = Cor, colour = Phenotype))+
    geom_freqpoly(bins = 20)+
    labs(x = "Correlation",
         y = "Frequency")
dev.off()


pdf(paste0(cor_dir,"Results/RNAseq_cor/Gene_cor_sig.pdf"), height = 3, width = 5)
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
dev.off()

pdf(paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_scatter.pdf"), height = 5, width = 7)
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

heat_mat <- heat_df %>% select(ACR, GFR) %>% as.matrix()
rownames(heat_mat) <- heat_df %>% pull(EnsID)

pdf(paste0(cor_dir,"Results/RNAseq_cor/15wk_RNAseq_cor_heatmap_sig.pdf"), height = 5, width = 3)
pheatmap(heat_mat, cutree_rows = 2, scale = "none",
  annotation_names_row = FALSE,
  show_rownames = FALSE)
dev.off()

hist(cor$acr_fdr)
