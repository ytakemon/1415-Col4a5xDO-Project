message("downloading devtools")
library(devtools)
message("downloading knitr")
library(knitr)
message("downloading DOQTL")
library(DOQTL)

setwd("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction")

message("loading grid files")
#load all grid files
grid.64K <- read.table("./reference.files/ref.genome_grid.64k.noYnoMT.txt", sep = "\t", header = T)
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

message("converting bp to Mb")
#change bp to Mb for 64k grid 
grid.64K[2] <- grid.64K[2] / 1000000

#create marker names for 64k grid
temp <- c("chr", "bp")
grid.64K$snp_name <- do.call(paste, c(grid.64K[temp], sep="_"))
grid.64K <- grid.64K[,c(4,1,2)]

message("loading civet files")
#load in sample list and file locations
civet.list <- read.table("./reconstruct/civet_list.txt", sep = "\t", header = F)
civet.files <- list.files("./reconstruct/gbrs.64k.interpol.reconstruction", full.names = T, pattern = ".tsv")

#R doesnt like it when things start with numbers, so chance to X1415.0000 format
civet.list <- as.character(civet.list[,1])
civet.list <- sub("^", "X", civet.list)
civet.list <- sub("-", ".", civet.list)

#construct individual dataframe with SNP names and proper founder names
founder.names <- c("A","B","C","D","E","F","G","H")

message("Adding snp names and CC founder names")
for (i in 1:length(civet.files)){
  tmp <- read.delim(civet.files[i], sep = "\t", header = T)
  rownames(tmp) <- grid.64K$snp_name 
  names(tmp) <- founder.names
  assign(civet.list[i], tmp)
}

message("Interpolating file to GM_snp")
#now all samples are in the correct format for interpolation 
for (i in 1:length(civet.list)){
  tmp <- eval(parse(text = civet.list[i])) #need to do this to be recognized as var.
  tmp2 <- interpolate.markers(tmp, grid.64K, GM_snps)
  write.table(tmp2, file = paste("./reconstruct/DOQTL.GM.interpol.reconstruction/", civet.list[i], ".GM.interpolated.genoprobs.tsv", sep = ""), sep = "\t")
}
