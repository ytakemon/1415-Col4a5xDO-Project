#source("http://bioconductor.org/biocLite.R")
#biocLite("DOQTL")
#library(devtools)
#install_github("dmgatti/DOQTL", force = TRUE)
library(DOQTL)

setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/Gillian/")
load("./GLB_possibleDOF1s.Rdata")


Gillian.probs <- probs
Gillian.names <- rownames(Gillian.probs)

for (i in 1:length(Gillian.names)){
  temp <- Gillian.probs[i,,]
  assign(Gillian.names[i], t(temp))
}

#extract markers from Gillian's data
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
GL_snps <- as.data.frame(rownames(GLB17))
rownames(GL_snps) <- GL_snps[,1]
GM_GL_snps <- intersect(rownames(GM_snps), rownames(GL_snps))
GM_GL_snps <- GM_snps[GM_GL_snps,]

#fix B6 distribution
#GL17 looks fine
Gillian.B6redist.probs <- Gillian.probs[-8,,]
for (i in 1:length(rownames(Gillian.B6redist.probs))){
  temp <- t(Gillian.B6redist.probs[i,,])
  temp[, 2] <- temp[, 2] - 0.5
  temp[,2][temp[,2] < 0] <- 0
  temp <- temp / rowSums(temp)
  Gillian.B6redist.probs[i,,] <- t(temp)
}
#add back G17
GLB17.array <- array(0, c(1, 8, 120789))
GLB17.array[1,,] <- Gillian.probs["GLB17",,]

Gillian.B6redist.probs <- abind(Gillian.B6redist.probs, GLB17.array, along = 1)
rownames(Gillian.B6redist.probs)[13] <- "GLB17"
#save B6 fixed data

#kinship matrix
K_GM_GL <- kinship.probs(Gillian.B6redist.probs, snps = GM_GL_snps, bychr= TRUE)


#create kinship map for each sample
load("./Gillian.B6redist.kinship.Rdata")
load("./Gillian.B6redist.probs.Rdata")

for (i in 1:length(rownames(Gillian.B6redist.probs))){
  tiff(paste0(rownames(Gillian.B6redist.probs)[i], ".reconst.tiff"))
  image(Gillian.B6redist.probs[i,,], main = rownames(Gillian.B6redist.probs)[i])
  dev.off()
  }

tiff("GillianDOF1like.tiff", width = 1000, height = 1000)
layout(matrix(1:13,3,5))
image(Gillian.B6redist.probs["GLB17",,], main = "GLB17")
image(Gillian.B6redist.probs["GLB173",,], main = "GLB173")
image(Gillian.B6redist.probs["GLB174",,], main = "GLB174")
image(Gillian.B6redist.probs["GLB175",,], main = "GLB175")
image(Gillian.B6redist.probs["GLB176",,], main = "GLB176")
image(Gillian.B6redist.probs["GLB177",,], main = "GLB177")
image(Gillian.B6redist.probs["GLB178",,], main = "GLB178")
image(Gillian.B6redist.probs["GLB179",,], main = "GLB179")
image(Gillian.B6redist.probs["GLB180",,], main = "GLB180")
image(Gillian.B6redist.probs["GLB181",,], main = "GLB181")
image(Gillian.B6redist.probs["GLB182",,], main = "GLB182")
image(Gillian.B6redist.probs["GLB183",,], main = "GLB183")
image(Gillian.B6redist.probs["GLB184",,], main = "GLB184")
dev.off()



################ Compare kinship of Gillian's 13 samples to ours----------------------------------------------------------------

#first compile X1415 reconstruction into 1 genoprob array
#sample_names <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/DOQTL.GM.interpol.reconstruction/X1415.names.only.txt", sep = "\t", header = F)
#X1415.0383 <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/DOQTL.GM.interpol.reconstruction/X1415.0383.GM.interpolated.genoprobs.tsv", sep = "\t", header = T)
#X1415.0170 <- read.delim("~/Desktop/Col4a5xDO/R_col4a5xdo/GM.interpol.genoprob/interpol.genoprob.0170.tsv", sep = "\t", header = T)
sample_names <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/DOQTL.GM.interpol.reconstruction/X1415.names.only.txt", sep = "\t", header = F)
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
load("./Gillian.data/Gillian_markers.Rdata") #GM_GL_snps
founder_names <- c("A","B","C","D","E","F","G","H")
snp_names <- rownames(GM_GL_snps)
all.tsv.filepath <- list.files("./B6.prob.redistribution/", full.names = T, pattern = ".tsv")


X1415.probs <- array(0, c(192, 8, 120789), dimnames = list(sample_names[,1], founder_names, snp_names))
for (i in 1:length(all.tsv.filepath)){
  temp <- read.delim(all.tsv.filepath[i], sep = "\t", header = T)
  temp <- temp[rownames(GM_GL_snps),] #subset to Gillian's marker #s
  X1415.probs[i,,] <- t(temp)
}
save(X1415.probs, file = "./Gillian.data/X1415.genoprob.Gillian.subet.Rdata")
#save file
#clear environment

#combine with Gillian's data
library(abind)
load("./Gillian.data/X1415.genoprob.Gillian.subet.Rdata")
load("./Gillian.data/Gillian_markers.Rdata")
load("./Gillian.data/Gillian.B6redist.probs.Rdata")


X1415.Gillian.combined.genoprobs <- abind(X1415.probs, Gillian.B6redist.probs, along = 1)
save(X1415.Gillian.combined.genoprobs, file = "./Gillian.data/X1415.Gillian.combiined.genoprobs.Rdata")

#Make kinship probability file
K_combined <- kinship.probs(X1415.Gillian.combined.genoprobs, snps = GM_GL_snps, bychr = T)
save( K_combined, file = "./Gillian.data/Combined_kinship.Rdata")

load("~/Desktop/Col4a5xDO/R_col4a5xdo/Gillian/Combined_kinship.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/Gillian/Gillian_markers.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/Gillian/X1415.Gillian.combiined.genoprobs.Rdata")
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/Gillian/")
#Map kinship heat map
#total
image(1:nrow(K_combined[[1]]), 1:ncol(K_combined[[1]]), K_combined[[1]][,ncol(K_combined[[1]]):1], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of 1415 and Gillian's data",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 20:0, las = 1)
#close up#1
image(160:205, 160:205, K_combined[[1]][160:205 , 205:160], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of 1415 and Gillian's data",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 20:0, las = 1)
#close up#2
image(25:50, 180:200, K_combined[[1]][25:50 , 200:180], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of 1415 and Gillian's data",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 20:0, las = 1)

#Isolating genotypes of possible duplicates:
#A1 160, 176
tiff("./A1.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[160,,], main = rownames(X1415.Gillian.combined.genoprobs)[160])
image(X1415.Gillian.combined.genoprobs[176,,], main = rownames(X1415.Gillian.combined.genoprobs)[176])   
dev.off()

#A2
tiff("./A2.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[161,,], main = rownames(X1415.Gillian.combined.genoprobs)[161])
image(X1415.Gillian.combined.genoprobs[177,,], main = rownames(X1415.Gillian.combined.genoprobs)[177])   
dev.off()     

#A3 & 4
tiff("./A3_4.tiff")
layout(matrix(1:3,1,3))
image(X1415.Gillian.combined.genoprobs[174,,], main = rownames(X1415.Gillian.combined.genoprobs)[174])
image(X1415.Gillian.combined.genoprobs[175,,], main = rownames(X1415.Gillian.combined.genoprobs)[175])  
image(X1415.Gillian.combined.genoprobs[199,,], main = rownames(X1415.Gillian.combined.genoprobs)[199])   
dev.off()

#A5
tiff("./A5.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[180,,], main = rownames(X1415.Gillian.combined.genoprobs)[180])
image(X1415.Gillian.combined.genoprobs[203,,], main = rownames(X1415.Gillian.combined.genoprobs)[203])   
dev.off()

#A6
tiff("./A6.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[181,,], main = rownames(X1415.Gillian.combined.genoprobs)[181])
image(X1415.Gillian.combined.genoprobs[204,,], main = rownames(X1415.Gillian.combined.genoprobs)[204])   
dev.off()

#A7
tiff("./A7.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[182,,], main = rownames(X1415.Gillian.combined.genoprobs)[182])
image(X1415.Gillian.combined.genoprobs[197,,], main = rownames(X1415.Gillian.combined.genoprobs)[197])   
dev.off()

#A8
tiff("./A8.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[183,,], main = rownames(X1415.Gillian.combined.genoprobs)[183])
image(X1415.Gillian.combined.genoprobs[200,,], main = rownames(X1415.Gillian.combined.genoprobs)[200])   
dev.off()

#A9
tiff("./A9.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[184,,], main = rownames(X1415.Gillian.combined.genoprobs)[184])
image(X1415.Gillian.combined.genoprobs[202,,], main = rownames(X1415.Gillian.combined.genoprobs)[202])   
dev.off()

#A10
tiff("./A10.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[185,,], main = rownames(X1415.Gillian.combined.genoprobs)[185])
image(X1415.Gillian.combined.genoprobs[201,,], main = rownames(X1415.Gillian.combined.genoprobs)[201])   
dev.off()

#A11
tiff("./A11.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[186,,], main = rownames(X1415.Gillian.combined.genoprobs)[186])
image(X1415.Gillian.combined.genoprobs[194,,], main = rownames(X1415.Gillian.combined.genoprobs)[194])   
dev.off()

#A12
tiff("./A12.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[187,,], main = rownames(X1415.Gillian.combined.genoprobs)[187])
image(X1415.Gillian.combined.genoprobs[193,,], main = rownames(X1415.Gillian.combined.genoprobs)[193])   
dev.off()

#A13
tiff("./A13.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[190,,], main = rownames(X1415.Gillian.combined.genoprobs)[190])
image(X1415.Gillian.combined.genoprobs[195,,], main = rownames(X1415.Gillian.combined.genoprobs)[195])   
dev.off()

#A14
tiff("./A14.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[191,,], main = rownames(X1415.Gillian.combined.genoprobs)[191])
image(X1415.Gillian.combined.genoprobs[196,,], main = rownames(X1415.Gillian.combined.genoprobs)[196])   
dev.off()

#A15
tiff("./A15.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[194,,], main = rownames(X1415.Gillian.combined.genoprobs)[194])
image(X1415.Gillian.combined.genoprobs[201,,], main = rownames(X1415.Gillian.combined.genoprobs)[201])   
dev.off()

#B1
tiff("./B1.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[33,,], main = rownames(X1415.Gillian.combined.genoprobs)[33])
image(X1415.Gillian.combined.genoprobs[189,,], main = rownames(X1415.Gillian.combined.genoprobs)[189])   
dev.off()

#B2
tiff("./B2.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[34,,], main = rownames(X1415.Gillian.combined.genoprobs)[34])
image(X1415.Gillian.combined.genoprobs[190,,], main = rownames(X1415.Gillian.combined.genoprobs)[190])   
dev.off()

#B3
tiff("./B3.tiff")
layout(matrix(1:2,1,2))
image(X1415.Gillian.combined.genoprobs[34,,], main = rownames(X1415.Gillian.combined.genoprobs)[34])
image(X1415.Gillian.combined.genoprobs[195,,], main = rownames(X1415.Gillian.combined.genoprobs)[195])   
dev.off()


####Run correlation on ones that were flagged to look similar:
#A34
A34 <- diag(cor(t(X1415.Gillian.combined.genoprobs[175,,]), t(X1415.Gillian.combined.genoprobs[199,,])))
A5 <- diag(cor(t(X1415.Gillian.combined.genoprobs[180,,]), t(X1415.Gillian.combined.genoprobs[203,,])))
A6 <- diag(cor(t(X1415.Gillian.combined.genoprobs[181,,]), t(X1415.Gillian.combined.genoprobs[204,,])))
A7 <- diag(cor(t(X1415.Gillian.combined.genoprobs[182,,]), t(X1415.Gillian.combined.genoprobs[197,,])))
A8 <- diag(cor(t(X1415.Gillian.combined.genoprobs[183,,]), t(X1415.Gillian.combined.genoprobs[200,,])))
A9 <- diag(cor(t(X1415.Gillian.combined.genoprobs[184,,]), t(X1415.Gillian.combined.genoprobs[202,,])))
A10 <- diag(cor(t(X1415.Gillian.combined.genoprobs[185,,]), t(X1415.Gillian.combined.genoprobs[201,,])))
A11 <- diag(cor(t(X1415.Gillian.combined.genoprobs[186,,]), t(X1415.Gillian.combined.genoprobs[194,,])))
A12 <- diag(cor(t(X1415.Gillian.combined.genoprobs[187,,]), t(X1415.Gillian.combined.genoprobs[193,,])))
A13 <- diag(cor(t(X1415.Gillian.combined.genoprobs[190,,]), t(X1415.Gillian.combined.genoprobs[195,,])))
A14 <- diag(cor(t(X1415.Gillian.combined.genoprobs[191,,]), t(X1415.Gillian.combined.genoprobs[196,,])))

X1415.Gillian.correlation <- matrix(0, nrow = 11, ncol = 8)
  X1415.Gillian.correlation[1,] <- A34
  X1415.Gillian.correlation[2,] <- A5
  X1415.Gillian.correlation[3,] <- A6
  X1415.Gillian.correlation[4,] <- A7
  X1415.Gillian.correlation[5,] <- A8
  X1415.Gillian.correlation[6,] <- A9
  X1415.Gillian.correlation[7,] <- A10
  X1415.Gillian.correlation[8,] <- A11
  X1415.Gillian.correlation[9,] <- A12
  X1415.Gillian.correlation[10,] <- A13
  X1415.Gillian.correlation[11,] <- A14

  
cor_samples <- c("X1415.1014_GLB179", "X1415.1029_GLB183", "X1415.1035_GLB184", "X1415.1036_GLB177", "X1415.1043_GLB180", "X1415.1051_GLB182", "X1415.1054_GLB181", "X1415.1056_GLB174", "X1415.1083_GLB173", "X1415.1095_GLB175", "X1415.1117_GLB176")
founder_names <- c("A","B","C","D","E","F","G","H")
colnames(X1415.Gillian.correlation) <- founder_names
rownames(X1415.Gillian.correlation) <- cor_samples
X1415.Gillian.correlation <- as.data.frame(X1415.Gillian.correlation)
X1415.Gillian.correlation$Avg <- rowSums(X1415.Gillian.correlation) / 8
write.table(X1415.Gillian.correlation, file = "./X1415.Gillian.correlations.txt")


#############correlate 11 MUGA data to RNA-seq----------------------------------------------------------------
#subset 11 X1415 samples from above and subset its MUGA array data
#then compare them to 192 RNA-seq dta

sample_names <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/DOQTL.GM.interpol.reconstruction/X1415.names.only.txt", sep = "\t", header = F)
founder_names <- c("A","B","C","D","E","F","G","H")
all.tsv.filepath <- list.files("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/B6.prob.redistribution", full.names = T, pattern = ".tsv")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
test <- read.delim("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.1122.B6.prob.redistributed.genoprob.tsv", sep = "\t", header = T)
snp_names <- rownames(test)

X1415.genoprobs <- array(0, c(192, 8, 143064), dimnames = list(sample_names[,1], founder_names, snp_names))
for (i in 1:length(all.tsv.filepath)){
  temp <- read.delim(all.tsv.filepath[i], sep = "\t", header = T)
  X1415.genoprobs[i,,] <- t(temp)
}
save(X1415.genoprobs, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Rdata")

#dims in MUGA vs reconstruction differ so subset reconstruction to fit MUGA data
load("/hpcdata/ytakemon/Col4a5xDO/megamuga/newcol4a5/col4a5_probs.Rdata")
load("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Rdata")

snp_names <- colnames(probs[1,,])
snp_names <- as.data.frame(snp_names)
subset.X1415.genoprobs <- X1415.probs <- array(0, c(192, 8, 120930), dimnames = list(sample_names[,1], founder_names, snp_names[,1]))
for (i in 1:192){
  temp <- t(X1415.genoprobs[i,,])
  temp <- temp[snp_names[,1],] #subset reconstruction data to Dan's genoprob dataset 
  subset.X1415.genoprobs[i,,] <- t(temp)
}
save(subset.X1415.genoprobs, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Dan.subset.Rdata")


########################## MAKE MUGE VS RNA-SEQ CORRELATION------------------------------------------------------------------
load("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Dan.subset.Rdata")
load("/hpcdata/ytakemon/Col4a5xDO/megamuga/newcol4a5/col4a5_probs.Rdata")

# Set the sample IDs for probs to all be 4 digits long, consistent with pheno data
colnames(probs) = sub("I", "", colnames(probs)) 
spl = strsplit(rownames(probs), "\\.")
one = sapply(spl, "[", 1)
two = sapply(spl, "[", 2)
wh = which(nchar(two) == 3)
two[wh] = paste0("0", two[wh])
rownames(probs) = paste(one, two, sep = ".")


RNA_seq <- subset.X1415.genoprobs
MUGA <- probs
founder_names <- c("A","B","C","D","E","F","G","H")
colnames(MUGA) <- founder_names
colnames(RNA_seq) <- founder_names

#X1415.1014 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1014",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1014.cor <- temp
save(X1415.1014.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1014.cor.Rdata")

#X1415.1029 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1029",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1029.cor <- temp
save(X1415.1029.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1029.cor.Rdata")

#X1415.1035 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1035",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1035.cor <- temp
save(X1415.1035.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1035.cor.Rdata")

#X1415.1036 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1036",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1036.cor <- temp
save(X1415.1036.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1036.cor.Rdata")

#X1415.1043 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1043",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1043.cor <- temp
save(X1415.1043.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1043.cor.Rdata")

#X1415.1051 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1051",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1051.cor <- temp
save(X1415.1051.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1051.cor.Rdata")

#X1415.1054 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1054",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1054.cor <- temp
save(X1415.1054.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1054.cor.Rdata")

#X1415.1056 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1056",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1056.cor <- temp
save(X1415.1056.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1056.cor.Rdata")

#X1415.1084 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1084",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1084.cor <- temp
save(X1415.1084.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1084.cor.Rdata")

#X1415.1095 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1095",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1095.cor <- temp
save(X1415.1095.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1095.cor.Rdata")

#X1415.1117 vs RNA-seq----------------------------------------
temp <- matrix(0, nrow = 192, ncol = 8)
rownames(temp) <- rownames(RNA_seq)
colnames(temp) <- founder_names
for (i in 1:192){
  temp2 <- diag(cor(t(MUGA["X1415.1117",,]), t(RNA_seq[i,,])))
  temp[i,] <- temp2
}
temp <- as.data.frame(temp)
temp$AvgCor <- rowSums(temp) / 8 
X1415.1117.cor <- temp
save(X1415.1117.cor, file = "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.correlation/X1415.1117.cor.Rdata")

#visualization MUGA vs RNA-seq-----------------------------------------------------------
library(ggplot2)

#1415.1014
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1014.cor.Rdata")
tempx <- rownames(X1415.1014.cor)
tempy <- X1415.1014.cor$AvgCor
ggplot(X1415.1014.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1014") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1014.tiff")

#1415.1029
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1029.cor.Rdata")
tempx <- rownames(X1415.1029.cor)
tempy <- X1415.1029.cor$AvgCor
ggplot(X1415.1014.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1029") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1029.tiff")

#1415.1035
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1035.cor.Rdata")
tempx <- rownames(X1415.1035.cor)
tempy <- X1415.1035.cor$AvgCor
ggplot(X1415.1035.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1035") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1035.tiff")

#1415.1036
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1036.cor.Rdata")
tempx <- rownames(X1415.1036.cor)
tempy <- X1415.1036.cor$AvgCor
ggplot(X1415.1036.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1036") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1036.tiff")

#1415.1043
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1043.cor.Rdata")
tempx <- rownames(X1415.1043.cor)
tempy <- X1415.1043.cor$AvgCor
ggplot(X1415.1043.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1043") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1043.tiff")

#1415.1051
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1051.cor.Rdata")
tempx <- rownames(X1415.1051.cor)
tempy <- X1415.1051.cor$AvgCor
ggplot(X1415.1051.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1051") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1051.tiff")

#1415.1054
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1054.cor.Rdata")
tempx <- rownames(X1415.1054.cor)
tempy <- X1415.1054.cor$AvgCor
ggplot(X1415.1054.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1054") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1054.tiff")

#1415.1056
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1056.cor.Rdata")
tempx <- rownames(X1415.1056.cor)
tempy <- X1415.1056.cor$AvgCor
ggplot(X1415.1056.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1056") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1056.tiff")

#1415.1084
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1084.cor.Rdata")
tempx <- rownames(X1415.1084.cor)
tempy <- X1415.1084.cor$AvgCor
ggplot(X1415.1084.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1084") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1084.tiff")

#1415.1095
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1095.cor.Rdata")
tempx <- rownames(X1415.1095.cor)
tempy <- X1415.1095.cor$AvgCor
ggplot(X1415.1095.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1095") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1095.tiff")

#1415.1117
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1117.cor.Rdata")
tempx <- rownames(X1415.1117.cor)
tempy <- X1415.1117.cor$AvgCor
ggplot(X1415.1117.cor, aes(tempx, tempy)) +
  geom_point (aes(y = AvgCor)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Correlation between RNA-seq reconstruction and X1415.1117") +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave(filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.correlation/X1415.1117.tiff")

############## Gillian vs RNA-seq reconstruction------------------------------------------------

#Gillian vs RNA-seq reconstruction
load("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/Gillian.B6redist.probs.Rdata")
load("/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/X1415.genoprob.Gillian.subet.Rdata")

#Variables: 
#Gillian.B6redist.probs
#X1415.probs

Gillian.sample.names <- rownames(Gillian.B6redist.probs)
founder_names <- c("A","B","C","D","E","F","G","H")
output <- "/hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/Gillian.data/GLB.correlation/"

for (i in 1:13){
  temp <- matrix(0, nrow = 192, ncol = 8)
  rownames(temp) <- rownames(X1415.probs)
  colnames(temp) <- founder_names
  for (x in 1:192){
    temp2 <- diag(cor(t(Gillian.B6redist.probs[i,,]), t(X1415.probs[x,,])))
    temp[x,] <- temp2
  }
  temp <- as.data.frame(temp)
  temp$AvgCor <- rowSums(temp) / 8
  assign(Gillian.sample.names[i], temp)
  save( list = paste(Gillian.sample.names[i]), file = paste(output, "/", Gillian.sample.names[i], ".cor.Rdata", sep = "") )
}
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB17.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB173.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB174.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB175.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB176.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB177.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB178.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB179.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB180.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB181.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB182.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB183.cor.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB184.cor.Rdata")

varname <- ls()

for (i in 1:13){
  Sample <- get(varname[i])
  tempx <- rownames(Sample)
  tempy <- Sample$AvgCor
    ggplot(Sample, aes(tempx, tempy)) +
      geom_point(aes(y = AvgCor)) +
      theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
      ggtitle(paste("Correlation between RNA-seq reconstruction and ", varname[i], sep = "")) +
      xlab("RNA-seq reconstructed samples") +
      ylab("Pearson Correlation")
    ggsave(filename = paste("~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/", varname[i], ".tiff", sep = "\t"))
}

### overlay points and labels-----------------

#GLB173
Sample <- GLB173

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 0.5 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB173")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB173_labeled.tiff")

#GLB174
Sample <- GLB174

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 0.5 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB174")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB174_labeled.tiff")

#GLB175
Sample <- GLB175

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB175")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB175_labeled.tiff")

#GLB176
Sample <- GLB176

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB176")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB176_labeled.tiff")

#GLB177
Sample <- GLB177

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB177")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB177_labeled.tiff")

#GLB178
Sample <- GLB178

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB178")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB178_labeled.tiff")

#GLB179
Sample <- GLB179

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB179")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB179_labeled.tiff")

#GLB180
Sample <- GLB180

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB180")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB180_labeled.tiff")


#GLB181
Sample <- GLB181

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB181")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB181_labeled.tiff")


#GLB182
Sample <- GLB182

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB182")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB182_labeled.tiff")

#GLB183
Sample <- GLB183

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB183")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB183_labeled.tiff")

#GLB184
Sample <- GLB184

max.name <- rownames(Sample)[which.max(Sample$AvgCor)]
max.value <- max(Sample$AvgCor)
MAX <- subset(Sample, AvgCor == max.value)
tempx <- rownames(Sample)
tempy <- Sample$AvgCor
tempmax <- rownames(MAX)
tempmay <- MAX$AvgCor

ggplot(Sample, aes(tempx, tempy)) +
  geom_point(data = Sample, aes(y = AvgCor)) +
  geom_point(data = MAX, aes(x = tempmax, y = tempmay, color = "red")) +
  geom_text( data = MAX, aes(x = tempmax, y = tempmay,label = max.name, hjust = 1 , vjust = 2)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between RNA-seq reconstruction and GLB184")) +
  xlab("RNA-seq reconstructed samples") +
  ylab("Pearson Correlation")
ggsave( filename = "~/Desktop/Col4a5xDO/R_col4a5xdo/GLB.correlation.RNA-seq/GLB184_labeled.tiff")



