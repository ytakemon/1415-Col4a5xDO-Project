#After going through both RNA-seq  file, GeneSeek files, and Gillian's files we have identified the best samples to use
#in this analyais.
#We will be taking 158 GS samples, 11 Gillian GS files, 1 re-analyzed files, and 22 RNA-seq reconstructed files
#for a total of 192 genoprob data array. 
#note: as of 09/12/16 Dan is stil working on the haplotype reconstructon of the 1 re-analyzed file so I will add this once he is finished.

#few things that need fixing: 
#1. X1415.1091 from GS is actually X1415.0243, and vice versa
#2. X1415.1016 from GS is actually X1415.0914, and vice versa


setwd("/hpcdata/ytakemon/Col4a5xDO")

GS157 <- read.delim("./Sample_list/Best.sample.list/GS157.txt", sep = "\t", header = FALSE)
Gillian11 <- read.delim("./Sample_list/Best.sample.list/Gillian11.txt", sep = "\t", header = FALSE)
RS24 <- read.delim("./Sample_list/Best.sample.list/RS24.txt", sep = "\t", header = FALSE)
ReGS1 <- "X1415.1039"
ReGS1 <- as.data.frame(ReGS1)

total.sample.names <- read.delim("./Sample_list/Best.sample.list/total.sample.list202.txt", sep = "\t", header = FALSE)

#total number of samples should be 202 until Dan finishes his reanalysis which should give us total of 203 samples

load("./GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Rdata") #RS data
RS.genoprobs <- X1415.genoprobs
load("./GBRS_reconstruction/reconstruct/Gillian.data/Gillian.B6redist.probs.Rdata") #Gillian data
Gillian.genoprobs <- Gillian.B6redist.probs
load("./megamuga/newcol4a5/col4a5_probs.Rdata") #GS data
GS.genoprobs <- probs

#rename 11 Gillian samples to 1415 IDs
dimnames(Gillian.genoprobs)[[1]][1] <- "X1415.1084" 
dimnames(Gillian.genoprobs)[[1]][2] <- "X1415.1056"
dimnames(Gillian.genoprobs)[[1]][3] <- "X1415.1095"
dimnames(Gillian.genoprobs)[[1]][4] <- "X1415.1117"
dimnames(Gillian.genoprobs)[[1]][5] <- "X1415.1036"
dimnames(Gillian.genoprobs)[[1]][7] <- "X1415.1014"
dimnames(Gillian.genoprobs)[[1]][8] <- "X1415.1043"
dimnames(Gillian.genoprobs)[[1]][9] <- "X1415.1054"
dimnames(Gillian.genoprobs)[[1]][10] <- "X1415.1051"
dimnames(Gillian.genoprobs)[[1]][11] <- "X1415.1029"
dimnames(Gillian.genoprobs)[[1]][12] <- "X1415.1035"

#fix GS.genoprobs names
spl = strsplit(rownames(GS.genoprobs), "\\.")
one = sapply(spl, "[", 1)
two = sapply(spl, "[", 2)
wh = which(nchar(two) == 3)
two[wh] = paste0("0", two[wh])
rownames(GS.genoprobs) = paste(one, two, sep = ".")

#sort genoprobs and name list
GS.genoprobs <- GS.genoprobs[ order(rownames(GS.genoprobs)),, ]
RS.genoprobs <- RS.genoprobs[ order(rownames(RS.genoprobs)),, ]
Gillian.genoprobs <- Gillian.genoprobs[ order(rownames(Gillian.genoprobs)),, ]

######################################################################################################################################
#Genoprobs 1: 157GS + 24RS + 11Gillian = 192 samples (best.genoprobs.192)
#Genoprobs 2: Only 157 GS genoprobs (genoprobs.GS.157)
#Genoprobs 3: 182GS (genoprobs.GS.182)

GS157 <- as.character(sort(GS157[,1]))
RS24 <- as.character(sort(RS24[,1]))
Gillian11 <- as.character(sort(Gillian11[,1]))

#Genoprobs 1
#subset  to use in "best.genoprobs.192"
GS.genoprobs.157 <- GS.genoprobs[GS157,,]
RS.genoprobs.24 <- RS.genoprobs[RS24,,]
Gillian.genoprobs.11 <- Gillian.genoprobs[Gillian11,,]
#Gillian's data set has the samllest dims. So subset all dato to her genoprob grid
RS.genoprobs.24 <- RS.genoprobs.24[,,dimnames(Gillian.genoprobs)[[3]]]
GS.genoprobs.157 <- GS.genoprobs.157[,,dimnames(Gillian.genoprobs)[[3]]]
#combine array using abind, combine by 3rd dimension
library(abind)
temp <- abind(GS.genoprobs.157, RS.genoprobs.24, along = 1)
best.genoprobs.192 <- abind(temp, Gillian.genoprobs.11, along = 1)
best.genoprobs.192 <- best.genoprobs.192[ order( rownames(best.genoprobs.192)),,]
save(best.genoprobs.192, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/best.genoprobs.192.Rdata")

#Genoprobs 2: 
#subset to use in "genoprobs.GS.157"
GS.genoprobs.157 <- GS.genoprobs[GS157,,]
GS.genoprobs.157 <- GS.genoprobs.157[,,dimnames(Gillian.genoprobs)[[3]]] #subset to smallest dims that's Gillians grid
GS.genoprobs.157 <- GS.genoprobs.157[ order( rownames(GS.genoprobs.157)),,] #order samples
save(GS.genoprobs.157, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GS.genoprobs.157.Rdata")

#Genoprobs 3: 
#subset to use in "genoprobs.168"
GS.genoprobs.182 <- GS.genoprobs
GS.genoprobs.182 <- GS.genoprobs.182[,,dimnames(Gillian.genoprobs)[[3]]]#subset to smallest dims that's Gillians grid
GS.genoprobs.182 <- GS.genoprobs.182[ order( rownames(GS.genoprobs.182)),,] # order samples
save(GS.genoprobs.182, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GS.genoprobs.182.Rdata")






