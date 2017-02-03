#R/3.3.1
library(reshape2)
library(ggplot2)
library(grid)

setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
pheno <- read.delim("./Phenotype/1415_master_pheno.txt", sep = "\t", header = TRUE)

#clean data
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192 samples
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
#keep only interested columns
pheno <- pheno[,c("MouseID", "Sex", "C2", "ACR6WK", "ACR10WK", "ACR15WK")]

pheno_F <- pheno[pheno$Sex == "F",]
pheno_M <- pheno[pheno$Sex == "M",]

#ACR Females
pheno_F <- pheno_F[,c("MouseID", "ACR6WK", "ACR10WK", "ACR15WK")]
ggdata<- melt(pheno_F)
names(ggdata) <- c("MouseID", "ACR_time", "Value")
ggplot_F <- ggplot(ggdata, aes(Value, ..count.., fill = ACR_time, colour = ACR_time)) +
	geom_density( alpha = 0.1)+
	scale_x_continuous("Albumin to creatinine ratio (mg/g)", limit = c(0 , 15000)) +
	scale_y_continuous(" ") +
	labs( title = "Females") +
	theme( plot.title = element_text(hjust = 0.5), 
			panel.background = element_rect(fill = "white", colour = "black"),
			axis.text.x = element_blank(), 
			axis.ticks.x = element_blank()
			) +
	coord_flip()

pheno_M <- pheno_M[,c("MouseID", "ACR6WK", "ACR10WK", "ACR15WK")]
ggdata<- melt(pheno_M)
names(ggdata) <- c("MouseID", "ACR_time", "Value")
ggplot_M <- ggplot(ggdata, aes(Value, ..count.., fill = ACR_time, colour = ACR_time)) +
	geom_density( alpha = 0.1)+
	scale_x_continuous("Albumin to creatinine ratio (mg/g)", limit = c(0 , 15000)) +
	scale_y_continuous(" ") +
	labs( title = "Males") +
	theme( plot.title = element_text(hjust = 0.5), 
			panel.background = element_rect(fill = "white", colour = "black"),
			axis.text.x = element_blank(), 
			axis.ticks.x = element_blank()
			) +
	coord_flip()

pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a5_fig2_ACR_GFR_dist_YUKA.pdf", width = 12, height = 6)
pushViewport(viewport(layout = grid.layout(1, 2)))
print(ggplot_F, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot_M, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

