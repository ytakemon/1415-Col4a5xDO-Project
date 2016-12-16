#create kinsihip to verify best.genoprob.192 sample are unique
library(DOQTL)
setwd("/hpcdata/ytakemon/Col4a5xDO")

load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/best.genoprobs.192.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata")) #GM_snps

K.probs  <- kinship.probs(best.genoprobs.192, snps = GM_snps, bychr = T)
save(K.probs, file = "./GBRS_reconstruction/reconstruct/best.compiled.genoprob/k.best.probs192.Rdata")

##############################################################################

#this part is done on the Rstudio
setwd("~/Desktop/Col4a5xDO/R_col4a5xdo/best.genoprob.192/")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/best.genoprob.192/k.best.probs192.Rdata")
load("~/Desktop/Col4a5xDO/R_col4a5xdo/best.genoprob.192/best.genoprobs.192.Rdata")

image(1:192, 1:192, K.probs[[1]][,1:192], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship map of best.genoprobs.192",
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 10 * 0:20, labels = 10 * 0:20, las = 1)
axis(side = 1, at = 10 * 0:20, labels = 10 * 0:20, las = 1)

#looks good
#ready for qtl analysis
