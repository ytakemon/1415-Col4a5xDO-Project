#Compiled Bayesian Intervals
#Bayesian intervals are derived from qtl maps 

#load GFR C2 model genoprobs data
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.GFR.log.C2.192.Rdata" ) # qtl.GFR.log.C2.192
qtl <- qtl.GFR.log.C2.192

#create interval for chr 2, 4, 7 ,9 ,10 ,15 , 19
interval = bayesint(qtl, chr = 2)
knitr::kable(interval)


#load Alb6wk genoprobs data
#create interval for chr: 2, 11, 12

#load Alb10wk genoprobs data
#create interval for chr: 2, 4, 11, 13
#load Alb15wk genoprobs data
#create interval for chr: 15