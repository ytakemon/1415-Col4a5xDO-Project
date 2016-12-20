#source and laod DOQTL package: this only needs to happen once
#source("http://bioconductor.org/biocLite.R")
#biocLite("DOQTL")
#library(devtools)
#library(knitr)
#install_github("dmgatti/DOQTL", force = TRUE)
library(DOQTL)
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/")

#load all grid file, Dan's 1415 genoprobs, GM_SNPs
grid_64k_noYnoMT <- read.table( "~/Desktop/Col4a5xDO/R_col4a5xdo/reconstruct.genoprob/ref.genome_grid.64k.noYnoMT.txt", sep = "\t", header = T)
temp <- c("chr", "bp")
grid_64k_noYnoMT$snp_name <- do.call(paste, c(grid_64k_noYnoMT[temp], sep="_"))


load("~/Desktop/Col4a5xDO/R_col4a5xdo/new_col4a5_probs.Rdata")

load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

#laod reconstructed samples
sample_names <- read.table("~/Desktop/Col4a5xDO/R_col4a5xdo/reconstruct.genoprob/reconstruct_list.txt", sep = "\t", header = F )
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/reconstruct.genoprob/")
list.reconst.XX.interp.files <- list.files("./noYnoMT.interpol/", full.names = T, pattern = ".tsv")

#combined.reconst.genoprobs <- lapply(list.reconst.XY.files, read.table, sep = "\t")

#make 3D array compatible with DOQTL package
#define dimnames
sample_names <- sample_names[[1]]
sample_names <- sub("^", "X", sample_names) #format Dan chose
sample_names <- sub("-", ".", sample_names) #format Dan chose
snp_names <- grid_64k_noYnoMT$snp_name
founder_names <- c("A","B","C","D","E","F","G","H")


#check dimentions of each reconstructed sample (just once to check)-------------------
reconst.count.XX <- array(0, c(33,2), dimnames = list(sample_names))
colnames(reconst.count.XX) <- c("Mrks","CC")
reconst.count.XX <-reconst.count.XX[order(row.names(reconst.count.XX)),]
for (i in 1:length(list.reconst.XX.interp.files)) {
  temp = read.delim(list.reconst.XX.interp.files[i], sep = "\t")
  temp2 = dim(temp)
  reconst.count.XX[i,] = temp2
}#dims are 64000 by 8

#subset out 44910 marker names
test <- read.delim(list.reconst.XX.interp.files[1], sep = "\t")

#Create a combined genoprobs data------------------------------------------------------
reconst.probs <- array(0, c(33,8,64000), dimnames = list(sample_names, founder_names, snp_names))
reconst.probs <- reconst.probs[order(row.names(reconst.probs)),,] #reorder 
for (i in 1:length(list.reconst.XX.interp.files)){
  temp = read.delim(list.reconst.XX.interp.files[i], sep = "\t")
  reconst.probs[i,,] <- t(temp)
}#dim of reconstructed genoprobs is 33, 8 , 64000

#create kinship probs
K_recon = kinship.probs(reconst.probs, snps = recon.snp, bychr = TRUE)

#Kinship mapping
image(1:nrow(K_recon[[1]]), 1:ncol(K_recon[[1]]), K_recon[[1]][,ncol(K_recon[[1]]):1], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship of reconstructed 1415 Col4a5xDO", 
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 20 * 0:7, labels = 20 * 7:0, las = 1)


#create kinship map for each sample---------------------------------------------------
for (i in 1:33){
    tiff(paste0(rownames(reconst.probs)[i], ".reconst.tiff"))
    image(reconst.probs[i,,], main = rownames(reconst.probs)[i])
    dev.off()
}



tiff("reconst.duplicate.samples.tiff")
layout(matrix(1:6,1,2))
image(reconst.probs["X1415.0176",,], main = "X1415.0176")
image(reconst.probs["X1415.0261",,], main = "X1415.0261")
dev.off()

#B6 like samples
tiff("reconst.B6like.samples.tiff")
layout(matrix(1:6,2,3))
image(reconst.probs["X1415.1014",,], main = "X1415.1014")
image(reconst.probs["X1415.1029",,], main = "X1415.1029")
image(reconst.probs["X1415.1035",,], main = "X1415.1035")
image(reconst.probs["X1415.1043",,], main = "X1415.1043")
image(reconst.probs["X1415.1051",,], main = "X1415.1051")
image(reconst.probs["X1415.1054",,], main = "X1415.1054")
dev.off()