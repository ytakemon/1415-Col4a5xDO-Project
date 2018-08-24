library(DOQTL)
library("biomaRt")
setwd("/hpcdata/ytakemon/Col4a5xDO")

# set up biomaRt
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

#	load files
load("./QTL/qtl.GFR.log.C2.192.Rdata")
load("./QTL/qtl.log.Alb6WK.192.Rdata")
load("./QTL/qtl.log.Alb10WK.192.Rdata")
load("./QTL/qtl.log.Alb15WK.192.Rdata")
#load("QTL/qtl_delta_ACR15_6.192.Rdata")
load("./QTL/perms.1000.qtl.GFR.log.C2.192.Rdata")
load("./QTL/perms.1000.qtl.log.Alb6WK.192.Rdata")
load("./QTL/perms.1000.qtl.log.Alb10WK.192.Rdata")
load("./QTL/perms.1000.qtl.log.Alb15WK.192.Rdata")
#load("./QTL/perms.1000.qtl.delta_ACR15_6.192.Rdata")

#	Create threshold by permutations
thr.1000.qtl.logAlb.6wk <- get.sig.thr( perms.1000.qtl.log.Alb6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.10wk <- get.sig.thr( perms.1000.qtl.log.Alb10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.15wk <- get.sig.thr( perms.1000.qtl.log.Alb15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
#thr.1000.qtl.deltaACR15_6 <- get.sig.thr(perms.1000.qtl.deltaACR15_6.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr  = FALSE)

## Threshold : C2
# Calcualte threshold
thr.1000.qtl.logGFR.6wk <- get.sig.thr( perms.1000.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
# show threshold
thr.1000.qtl.logGFR.6wk
#    0.05      0.1     0.63
#7.340369 6.945842 4.659221

# Generate gene list in invterval:
bayes <- bayesint(qtl, 19, 0.95)
interval <- bayes$pos * 1e6
getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "coding"))



# Plot QTL: GFR
#	remove X chr from qtl object
qtl <- qtl.GFR.log.C2.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Figure4.1.qtl.log.C2.GFR.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.C2.GFR QTL perms.1000.noX")
dev.off()





## Alb6wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb6WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Figure4.2.qtl.log.Alb6wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.6WK QTL perms.1000.noX")
dev.off()

## Alb10wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb10WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Figure4.3.qtl.log.Alb10wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.10wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.10WK QTL perms.1000.noX")
dev.off()

## Alb10wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb15WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Figure4.4.qtl.log.Alb15wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.15wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.15WK QTL perms.1000.noX")
dev.off()

## delta ACR15_6wk QLT plot
# remove X chr from qtl object
qtl <- qtl.deltaACR15_6.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Figure4.5.qtl.delta.ACR15_6.noX.pdf", width = 12, height = 6)
plot1 <- plot(qtl, sig.thr = thr.1000.qtl.deltaACR15_6, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO delta.ACR15_6 QTL perms.1000.noX")
dev.off()

## Plot all qtl into one pdf document
pushViewport(viewport( layout = grid.layout(5,1)))
viewport(layout.pos.row = 1, layout.pos.col = 1))




pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/PCA_M_pc12345.pdf", width = 6, height = 6)
pushViewport(viewport( layout = grid.layout(2,2)))
print(ggplot12, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot32, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ggplot34, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(ggplot54, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

library(grid)
library(gridBase)
plot.new()
vp1 <- viewport(x=0,y=0.5,width=0.5, height=0.5, just = c("left", "bottom"))
vp2 <- viewport(x=0.5,y=0,width=0.5, height=0.5, just = c("left", "bottom"))
pushViewport(vp1)
grid.rect()
grid.text("vp1", 0.5, 0.5)
par(new=TRUE, fig=gridFIG())
plot(1,2)
upViewport()
pushViewport(vp2)
grid.rect()
grid.text("vp2", 0.5, 0.5)
par(new=TRUE, fig=gridFIG())
plot(1,2)
