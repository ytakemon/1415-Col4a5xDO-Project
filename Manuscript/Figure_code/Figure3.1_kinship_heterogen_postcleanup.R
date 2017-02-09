#R/3.3.1
library(DOQTL)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

setwd("/hpcdata/ytakemon/Col4a5xDO")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/GM_snps.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/kinship/k.best.probs192.Rdata")

probs <- K.probs[[1]]
probs <- melt(probs)
names(probs)[3] <- "Correlation"
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Figure3.1_kinship_heterogen_postcleanup.pdf", width = 7, height = 7)
	qplot(x=Var1, y=Var2, data=probs, fill=Correlation, geom="tile") +
	scale_x_discrete("Col4a5 x D0 F1 samples") +
	scale_y_discrete("Col4a5 x D0 F1 samples") +
	theme(	axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) +
   scale_fill_gradient(limits=c(0, 1), low ="red", high ="white")
dev.off()

