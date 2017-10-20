library(dplyr)
setwd("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob")
load("genoprobs/best.genoprobs.192.Rdata")
load("qtl/RNA/Complete_eQTL_Rdata/ENSMUSG00000076147.Mir694.eQTL.rankZ.tpm.Rdata")
load("GM_snps.Rdata")
load("sex.covar.Rdata")
load("RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load("EnsemblID_GRCm38.p4.Rdata")
pheno <- read.csv("pheno.csv")

Ensembl_mouseID <- Ensembl_mouseID[Ensembl_mouseID$chromosome_name %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, "X"),]
Ensembl_mouseID <- arrange(Ensembl_mouseID, ensembl_gene_id) %>% select(ensembl_gene_id, mgi_symbol, chromosome_name, start_position, end_position)
