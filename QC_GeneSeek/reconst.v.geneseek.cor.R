#R/3.3.1
#To compare all RNA-seq reconstructed genoprobs to GeneSeek Muga genoprobs array
#To identify incongurencies between samples before beginning analysis and get our refunds for their mistake

setwd("/hpcdata/ytakemon/Col4a5xDO")


##############################
#Compile b6.prob genoprobs into one array

sample_names <- read.delim("./GBRS_reconstruction/reconstruct/DOQTL.GM.interpol.reconstruction/X1415.names.only.txt", sep = "\t", header = F)
founder_names <- c("A","B","C","D","E","F","G","H")
all.tsv.filepath <- list.files("./GBRS_reconstruction/reconstruct/B6.prob.redistribution", full.names = T, pattern = ".tsv")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
temp <- read.delim("./GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.1122.B6.prob.redistributed.genoprob.tsv", sep = "\t", header = T)
snp_names <- rownames(temp)

X1415.genoprobs <- array(0, c(192, 8, 143064), dimnames = list(sample_names[,1], founder_names, snp_names))
for (i in 1:length(all.tsv.filepath)){
  temp <- read.delim(all.tsv.filepath[i], sep = "\t", header = T)
  X1415.genoprobs[i,,] <- t(temp)
}
save(X1415.genoprobs, file = "./GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Rdata")


##############################
#Createt kinship map for RNA-seq reconstructed data

library(DOQTL)
library(abind)
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata")) #GM_snps
load("./GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Rdata") #X1415.genoprobs

K_X1415.RNA <- kinship.probs(X1415.genoprobs, snps = GM_snps, bychr = T)
save(K_X1415.RNA, file = "./GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.kinship.Rdata")

#map on desktop Rstudio
load("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.B6.redist/X1415.kinship.Rdata")
image(1:192, 1:192, K_X1415.RNA[[1]][,1:192], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of RNA reconstruct post B6.prob",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
axis(side = 1, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
#close up 1:
image(189:192, 30:40, K_X1415.RNA[[1]][189:192, 30:40], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of RNA reconstruct post B6.prob",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
axis(side = 1, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
#1: 189,33
#2: 190:34
#close up 2:
image(170:180, 155:165, K_X1415.RNA[[1]][170:180, 155:165], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of RNA reconstruct post B6.prob",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
axis(side = 1, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
#3: 176,160
#4: 177:161
#map pairs of samples
#1: 189,33
tiff("./pair1.comarison.tiff")
layout(matrix(1:2,1,2))
image(X1415.B6.Dan.subset.genoprob[189,,], main = rownames(X1415.B6.Dan.subset.genoprob)[189])
image(X1415.B6.Dan.subset.genoprob[33,,], main = rownames(X1415.B6.Dan.subset.genoprob)[33])
dev.off()

#2: 190:34
tiff("./pair2.comarison.tiff")
layout(matrix(1:2,1,2))
image(X1415.B6.Dan.subset.genoprob[190,,], main = rownames(X1415.B6.Dan.subset.genoprob)[190])
image(X1415.B6.Dan.subset.genoprob[34,,], main = rownames(X1415.B6.Dan.subset.genoprob)[34])
dev.off()

#3: 176,160
tiff("./pair3.comarison.tiff")
layout(matrix(1:2,1,2))
image(X1415.B6.Dan.subset.genoprob[176,,], main = rownames(X1415.B6.Dan.subset.genoprob)[176])
image(X1415.B6.Dan.subset.genoprob[160,,], main = rownames(X1415.B6.Dan.subset.genoprob)[160])
dev.off()

#4: 177:161
tiff("./pair4.comarison.tiff")
layout(matrix(1:2,1,2))
image(X1415.B6.Dan.subset.genoprob[177,,], main = rownames(X1415.B6.Dan.subset.genoprob)[177])
image(X1415.B6.Dan.subset.genoprob[161,,], main = rownames(X1415.B6.Dan.subset.genoprob)[161])
dev.off()

#everything looks good.


##############################
#Take each MUGA data and compare to 192 RNA-seq reconstructed samples

#dims in MUGA vs reconstruction differ so subset reconstruction to fit MUGA data
load("./megamuga/newcol4a5/col4a5_probs.Rdata")
load("./GBRS_reconstruction/reconstruct/B6.prob.redistribution/X1415.genoprob.Rdata")

X1415.B6.Dan.subset.genoprob <- X1415.genoprobs[,,dimnames(probs)[[3]]]
save(X1415.B6.Dan.subset.genoprob, file = "./GBRS_reconstruction/reconstruct/B6.prob.redistribution/B6.Dan.genoprob.subset.Rdata")

#Compare:
load("./megamuga/newcol4a5/col4a5_probs.Rdata")
load("./GBRS_reconstruction/reconstruct/B6.prob.redistribution/B6.Dan.genoprob.subset.Rdata")
Muga <- probs
RNA <- X1415.B6.Dan.subset.genoprob

#need to fix nameing of Muga data
founder_names <- c("A","B","C","D","E","F","G","H")
colnames(Muga) <- founder_names 

spl = strsplit(rownames(Muga), "\\.")
one = sapply(spl, "[", 1)
two = sapply(spl, "[", 2)
wh = which(nchar(two) == 3)
two[wh] = paste0("0", two[wh])
rownames(Muga) = paste(one, two, sep = ".")
Muga <- Muga[ order(rownames(Muga)),, ]

#create table of correlations comparing MUGA to RNA-reconstruction data
Muga_names <- rownames(Muga)
RNA_names <- rownames(RNA)
output <- "./GBRS_reconstruction/reconstruct/Muga.RNA.cor/"

temp <- matrix(0, nrow = 182, ncol = 192)
rownames(temp) <- Muga_names
colnames(temp) <- RNA_names

for (i in 1:182){
  for (x in 1:192){
    temp[i,x] <- cor(as.vector(Muga[i,,]), as.vector(RNA[x,,]))
  }
}
All.Muga.AvgCor <- as.data.frame(temp)

save(All.Muga.AvgCor, file = paste(output, "All.Muga.AvgCor.Rdata", sep = ""))
write.table(All.Muga.AvgCor, file = "./GBRS_reconstruction/reconstruct/Muga.RNA.cor/All.Muga.RNA.cor.txt", sep = "\t")

#make new col with Max covar
load("./GBRS_reconstruction/reconstruct/Muga.RNA.cor/All.Muga.AvgCor.Rdata")
All.Muga.AvgCor <- as.data.frame(t(All.Muga.AvgCor))
max <- apply(All.Muga.AvgCor, 1, max)
max <- as.data.frame(max)
write.table(max, file = "./GBRS_reconstruction/reconstruct/Muga.RNA.cor/MaxCor.txt", sep ="\t")
All.Muga.AvgCor$MaxCor <- max[,1]

##############################

#map all correlations for each of the 192 RNA-seq data
library(ggplot2)
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/X1415.B6.redist/")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/1415Muga.RNA.cor/All.Muga.AvgCor.Rdata")

for (i in 1:192){
sample.name <- colnames(All.Muga.AvgCor)[i]
tempx <- rownames(All.Muga.AvgCor)
tempy <- All.Muga.AvgCor[,i]

max.name <- rownames(All.Muga.AvgCor)[which.max(All.Muga.AvgCor[,i])]
max.value <- max(All.Muga.AvgCor[,i])
max.data <- subset(All.Muga.AvgCor, tempy == max.value)
max.x <- rownames(max.data)
max.y <- max.data[,i]


ggplot(All.Muga.AvgCor, aes(tempx, tempy)) +
  geom_point( data = All.Muga.AvgCor, aes( y = tempy)) +
  geom_point( data = max.data, aes(x = max.x, y = max.y, color = "red")) +
  geom_text( data = max.data, aes(x = max.x, y = max.y, label = max.name, hjust = 0.5, vjust = 2 )) +
  theme( axis.text.x = element_text( angle = -90, hjust = 0)) +
  ggtitle(paste("Correlation between reconstructed ", sample.name, " and GeneSeek MUGA samples")) +
           xlab("GeneSeek MUGA samples") +
           ylab("Pearson Correlation")
ggsave( filename = paste("~/Desktop/Col4a5xDO/R_col4a5xdo/1415Muga.RNA.cor/", sample.name, ".GS.correlation.tiff", sep = ""))
}


