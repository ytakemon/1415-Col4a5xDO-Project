################################################################################################################################
## Identifying trans-eQTL for Formin 1 (Fmn1)
## Yuka Takemon
## Created: 12/14/16
## Last modified: 12/19/16

## This script will identify trans-eQTL for Formin 1 by using the eQTL scanone objects created by eQTL_all.R Each scanone eQTL 
##	object will be assessed for significant peaks at Fmn1, more specifically at  "JAX00098823" or "rs27466023". Sigificance will
## 	be determined by their LOD scores above the pre-calulated LOD scores:
## > threshold
## sig.level         A        X
## 0.05 		8.104061 9.206618
## 0.1  		7.632182 9.048455
## 0.63 		6.047843 7.149671
## since we have 3 levels of significance we will be make 3 lists one for each threshold to better narrow down candidates that 
## have a trans-eQTL at Fmn1 ( "JAX00098823" or "rs27466023")

## Note:
## This script is tailored specifically for FMN1 however by changing marker parameters, we can apply this to find trans-eQTL for
## any genes. If I find time later, I can alter the script to make it more adaptable.

#	load library
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

#	Helperfuntions used:
source("substrR_L.R")

#	Input data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_rankZ_tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/EnsemblID_GRCm38.p4.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.RNA.rankZ.tpm.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/Threshold.perms.1000.qtl.RNA.rankZ.tpm.Rdata")
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

#Create empty array
Marker <- "JAX00098823" #define marker of interest
Fmn1_trans_eQTL <- array(0, c(length(All_eqtl_files), 2), dimnames = list (Ensembl_mouseID$ensembl_gene_id, c("MGI_ID", "LOD_score_at_JAX00098823")))
Fmn1_trans_eQTL <- as.data.frame(Fmn1_trans_eQTL)
Fmn1_trans_eQTL$MGI_ID <- Ensembl_mouseID$mgi_symbol

for (i in 1:length(All_eqtl_files)){
	#	load saved qtl object and assing to new variable name, temp, for easier manipulation
	rm(list = ls(pattern = "ENSMUSG"))	# reset environment clear of objects starting with values "ENSMUSG"
	load(All_eqtl_files_dir[i])
	temp <- get(ls(pattern = "ENSMUSG")) 
	#	temp$lod$A is a data frame with colnames: marker & lod	
	Fmn1_trans_eQTL$LOD_score_at_JAX00098823[i] <- temp$lod$A[Marker, "lod"]
	print(i)
	print(Fmn1_trans_eQTL$LOD_score_at_JAX00098823[i])
}

save(Fmn1_trans_eQTL, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Fmn1_trans_eQTL.Rdata")
#subset Fmn1_trans_eQTL by significance level  (something wrong with rownames of Fmn1)

Fmn1_trans_eQTL_Sig.0.05 <- Fmn1_trans_eQTL[Fmn1_trans_eQTL$LOD_score_at_JAX00098823 > 8.104061,]

Fmn1_trans_eQTL_Sig.0.1 <- Fmn1_trans_eQTL[Fmn1_trans_eQTL$LOD_score_at_JAX00098823 > 7.632182,]

Fmn1_trans_eQTL_Sig.0.63 <- Fmn1_trans_eQTL[Fmn1_trans_eQTL$LOD_score_at_JAX00098823 > 6.047843,]

write.table(Fmn1_trans_eQTL, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Fmn1_trans_eQTL_all.txt", sep="\t")
write.table(Fmn1_trans_eQTL_Sig.0.05, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Fmn1_trans_eQTL_Sig.0.05.txt", sep="\t")
write.table(Fmn1_trans_eQTL_Sig.0.1, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Fmn1_trans_eQTL_Sig.0.1.txt", sep="\t")
write.table(Fmn1_trans_eQTL_Sig.0.63, "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/Fmn1_trans_eQTL_Sig.0.63.txt", sep="\t")


