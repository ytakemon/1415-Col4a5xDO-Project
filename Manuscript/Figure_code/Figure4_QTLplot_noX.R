library(DOQTL)
library(ggplot2)
library(reshape2)
library(knitr)
setwd("/hpcdata/ytakemon/Col4a5xDO")
#	load files
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.GFR.log.C2.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15WK.192.Rdata")
#	Create threshold by permutations
thr.1000.qtl.logGFR.6wk <- get.sig.thr( perms.1000.qtl.GFR.log.C2.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.6wk <- get.sig.thr( perms.1000.qtl.log.Alb6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.10wk <- get.sig.thr( perms.1000.qtl.log.Alb10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.15wk <- get.sig.thr( perms.1000.qtl.log.Alb15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

## Alb6wk QTL plot
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