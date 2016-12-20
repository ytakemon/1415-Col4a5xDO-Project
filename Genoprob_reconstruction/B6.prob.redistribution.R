setwd("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction")
interpolated.genoprobs.dir <- "./reconstruct/DOQTL.GM.interpol.reconstruction/"

print(commandArgs(trailingOnly = TRUE))
sample.name <- commandArgs(trailingOnly = TRUE)
i <- read.delim( paste( interpolated.genoprobs.dir, sample.name, ".GM.interpolated.genoprobs.tsv", sep = ""), sep = "\t", header = T) 

#redistribute B6 probability 
i[,"B"] <- i[,"B"] - 0.5
i$B[i$B < 0] <- 0 

i <- i / rowSums(i)

range(rowSums(i))

#save data frame
write.table(i, file = paste("./reconstruct/B6.prob.redistribution/", sample.name, ".B6.prob.redistributed.genoprob.tsv", sep = ""), sep = "\t")

