library(pacman)
p_load(tidyverse, broom)

cor_dir <- "/projects/marralab/ytakemon_prj/Col4a5/"
load(paste0(cor_dir,"Data/RNAseq/RNA_seq_tpm.Rdata"))
load(paste0(cor_dir,"Data/Col4a5xDO_192data_YT.Rdata"))

RNA_seq <- RNA_seq %>% as_tibble(., rownames = "MouseID")
write_csv(RNA_seq, "~/Desktop/GBRS_processed_tpm.csv")

sample_annot <- Pheno %>% select(MouseID, Sex, SireGeneration)
write_csv(sample_annot, "~/Desktop/sample_annot.csv")
