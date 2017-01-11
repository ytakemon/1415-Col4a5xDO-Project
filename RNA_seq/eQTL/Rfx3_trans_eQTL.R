################################################################################################################################
## Identifying trans-eQTL for Rfx3
## Yuka Takemon
## Created: 01/10/17

## OBJECTIVE
## This script will identify trans-eQTL for Rfx3 by using the eQTL scanone objects created by eQTL_all.R Each scanone eQTL 
##	object will be assessed for significant peaks at Rfx3, more specifically at  "JAX00098823" or "rs27466023". Sigificance will
## 	be determined by their LOD scores above the pre-calulated LOD scores:
## > threshold
## sig.level         A        X
## 0.05 		8.104061 9.206618
## 0.1  		7.632182 9.048455
## 0.63 		6.047843 7.149671
## since we have 3 levels of significance we will be make 3 lists one for each threshold to better narrow down candidates that 
## have a trans-eQTL at Rfx3

## sessionInfo()
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## locale:
## [1] C
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets
## [8] methods   base
## other attached packages:
## [1] RColorBrewer_1.1-2   knitr_1.11           DOQTL_1.0.0
## [4] AnnotationDbi_1.28.2 GenomeInfoDb_1.2.5   IRanges_2.0.1
## [7] S4Vectors_0.4.0      Biobase_2.26.0       BiocGenerics_0.12.1
## [10] RSQLite_1.0.0        DBI_0.3.1
## loaded via a namespace (and not attached):
##  [1] Biostrings_2.34.1      GenomicRanges_1.18.4   MUGAExampleData_1.0.0
##  [4] QTLRel_0.2-14          RCurl_1.95-4.7         RUnit_0.4.30
##  [7] Rsamtools_1.18.3       XML_3.98-1.3           XVector_0.6.0
## [10] annotate_1.44.0        annotationTools_1.40.0 biomaRt_2.22.0
## [13] bitops_1.0-6           corpcor_1.6.8          gdata_2.17.0
## [16] gtools_3.5.0           hwriter_1.3.2          mclust_5.1
## [19] org.Hs.eg.db_3.0.0     org.Mm.eg.db_3.0.0     tools_3.1.1
## [22] xtable_1.8-0           zlibbioc_1.12.0

## List of data saved from this script (time savers for reanalysis)
## Rdata:
## save(Rfx3_trans_eQTL, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL.Rdata")
## Tables:
## write.table(Rfx3_trans_eQTL, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_all.txt", sep="\t")
## write.table(Rfx3_trans_eQTL_Sig.0.05, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_Sig.0.05.txt", sep="\t")
## write.table(Rfx3_trans_eQTL_Sig.0.1, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_Sig.0.1.txt", sep="\t")
## write.table(Rfx3_trans_eQTL_Sig.0.63, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_Sig.0.63.txt", sep="\t")
## Plots:
################################################################################################################################

#module load R/3.1.1

#	load library
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#	Helperfuntions used:
#source("substrR_L.R")

#	Input data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")

#	Get all eQTL file paths
All_eqtl_files_dir <- list.files("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata", full.names = T)
All_eqtl_files <- sub("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/","",All_eqtl_files_dir)
All_eqtl_files <- substrLeft(All_eqtl_files, 18)
#missing <- setdiff(Ensembl_mouseID$ensembl_gene_id, All_eqtl_files), data missing because I removed all QTL that had NA due to depredicated 
#ensembl identifiers.
#on cadillac rm *.NA.* and rm *..*

#	formatting
RNA_seqZ <- as.data.frame(RNA_seqZ)
Ensembl_mouseID <- Ensembl_mouseID[!duplicated(Ensembl_mouseID$ensembl_gene_id),]
rownames(Ensembl_mouseID) <- make.names(Ensembl_mouseID[,1]) 
Ensembl_mouseID <- Ensembl_mouseID[ order( rownames(Ensembl_mouseID)),]
Ensembl_mouseID <- Ensembl_mouseID[ intersect(rownames(Ensembl_mouseID), All_eqtl_files), ]

#	Quick directories
qtl_dir <- "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/RNA/Complete_eQTL_Rdata/"
plot_dir <- "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/RNA_qtl/Complete_eQTL_plots/"

#	Check if snps exists
GM_snps["JAX00474835",] #marker identified with highest lod score in region
# looks up rsID  rs30867551 to make sure that it falls within Rfx3
# No find another 
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata")
GFR_qtl <- qtl.GFR.log.C2.192$lod$A
GFR_qtl <- GFR_qtl[GFR_qtl[,"chr"] == 19, ]
GFR_qtl <- GFR_qtl[GFR_qtl[,"lod"] > 6, ]
# Now try UNCHS047474
GM_snps["UNCHS047474",] #Looks good!

#Create empty array
Marker <- "UNCHS047474" #define marker of interest
Rfx3_trans_eQTL <- array(0, c(length(All_eqtl_files), 2), dimnames = list (Ensembl_mouseID$ensembl_gene_id, c("MGI_ID", "LOD_score_at_UNCHS047474")))
Rfx3_trans_eQTL <- as.data.frame(Rfx3_trans_eQTL)
Rfx3_trans_eQTL$MGI_ID <- Ensembl_mouseID$mgi_symbol

for (i in 1:length(All_eqtl_files)){
	#	load saved qtl object and assing to new variable name, temp, for easier manipulation
	rm(list = ls(pattern = "ENSMUSG"))	# reset environment clear of objects starting with values "ENSMUSG"
	load(All_eqtl_files_dir[i])
	temp <- get(ls(pattern = "ENSMUSG")) 
	#	temp$lod$A is a data frame with colnames: marker & lod	
	Rfx3_trans_eQTL$LOD_score_at_UNCHS047474[i] <- temp$lod$A[Marker, "lod"]
	print(i)
	print(Rfx3_trans_eQTL$LOD_score_at_UNCHS047474[i])
}

save(Rfx3_trans_eQTL, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL.Rdata")
#subset Rfx3_trans_eQTL by significance level  (something wrong with rownames of Rfx3)

Rfx3_trans_eQTL_Sig.0.05 <- Rfx3_trans_eQTL[Rfx3_trans_eQTL$LOD_score_at_UNCHS047474 > 8.104061,] #sig level 0.05 for A

Rfx3_trans_eQTL_Sig.0.1 <- Rfx3_trans_eQTL[Rfx3_trans_eQTL$LOD_score_at_UNCHS047474 > 7.632182,] #sig level 0.1 for A

Rfx3_trans_eQTL_Sig.0.63 <- Rfx3_trans_eQTL[Rfx3_trans_eQTL$LOD_score_at_UNCHS047474 > 6.047843,] #sig level 0.63 for A

write.table(Rfx3_trans_eQTL, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_all.txt", sep="\t")
write.table(Rfx3_trans_eQTL_Sig.0.05, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_Sig.0.05.txt", sep="\t")
write.table(Rfx3_trans_eQTL_Sig.0.1, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_Sig.0.1.txt", sep="\t")
write.table(Rfx3_trans_eQTL_Sig.0.63, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/eQTL/Rfx3_eQTL/Rfx3_trans_eQTL_Sig.0.63.txt", sep="\t")


