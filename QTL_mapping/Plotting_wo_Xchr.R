###############################################################################################################################
## Additional QTL plotting:
#	Plotting QTL as .pdf without X chromosome
#	Plotting allele effect of baysian interval
#	Plotting average of founder effect by founder
## Yuka Takemon
## 12/05/16

## Setup
library(DOQTL)
library(ggplot2)
library(reshape2)
library(knitr)
setwd("/hpcdata/ytakemon/Col4a5xDO")
#	load files
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb6WK.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb10WK.192.Rdata")
load ("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/qtl/qtl.log.Alb15WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb6WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb10WK.192.Rdata")
load("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/perm1000/perms.1000.qtl.log.Alb15WK.192.Rdata")
#	Create threshold by permutations
thr.1000.qtl.logAlb.6wk <- get.sig.thr( perms.1000.qtl.log.Alb6WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.10wk <- get.sig.thr( perms.1000.qtl.log.Alb10WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)
thr.1000.qtl.logAlb.15wk <- get.sig.thr( perms.1000.qtl.log.Alb15WK.192[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

#	If you don't know you bayesian interval
#interval = bayesint(qtl, chr = 2)
#knitr::kable(interval)

## Helper functions (http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions)
#	Summarizes data, handling within-subjects variables by removing inter-subject variability.
#	It will still work if there are no within-S variables.
#	Gives count, un-normed mean, normed mean (with same between-group mean),
#	standard deviation, standard error of the mean, and confidence interval.
#	If there are within-subject variables, calculate adjusted values using method from Morey (2008).
#	data: a data frame.
#	measurevar: the name of a column that contains the variable to be summariezed
#	betweenvars: a vector containing names of columns that are between-subjects variables
#	withinvars: a vector containing names of columns that are within-subjects variables
#	idvar: the name of a column that identifies each subject (or matched subjects)
#	na.rm: a boolean that indicates whether to ignore NA's
#	conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
    FUN=is.factor, FUN.VALUE=logical(1))
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  # Apply correction from Morey (2008) to the standard error and confidence interval  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                           FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
#	Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#	data: a data frame.
#	measurevar: the name of a column that contains the variable to be summariezed
#	groupvars: a vector containing names of columns that contain grouping variables
#	na.rm: a boolean that indicates whether to ignore NA's
#	conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
#	Norms the data within specified groups in a data frame; it normalizes each
#	subject (identified by idvar) so that they have the same mean, within each group
#	specified by betweenvars.
#	data: a data frame.
#	idvar: the name of a column that identifies each subject (or matched subjects)
#	measurevar: the name of a column that contains the variable to be summariezed
#	betweenvars: a vector containing names of columns that are between-subjects variables
#	na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    library(plyr)

    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
     .fun = function(xx, col, na.rm) {
        c(subjMean = mean(xx[,col], na.rm=na.rm))
      },
      measurevar,
      na.rm
    )

    # Put the subject means with original data
    data <- merge(data, data.subjMean)

    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                               mean(data[,measurevar], na.rm=na.rm)

    # Remove this subject mean column
    data$subjMean <- NULL

    return(data)
}
## Alb6wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb6WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb6wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.6wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.6WK QTL perms.1000.noX")
dev.off()

## Alb10wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb10WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb10wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.10wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.10WK QTL perms.1000.noX")
dev.off()

## Alb10wk QTL plot
#	remove X chr from qtl object
qtl <- qtl.log.Alb15WK.192
qtl$lod$X <- NULL
qtl$coef$X <- NULL
#	plot qtl
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/qtl.log.Alb15wk.noX.pdf", width = 10.0, height = 7.5)
plot(qtl, sig.thr = thr.1000.qtl.logAlb.15wk, sig.col = c("red", "orange", "chartreuse"), main = "Col4a5xDO log.Alb.15WK QTL perms.1000.noX")
dev.off()

## Plotting narrower allele effect Alb 10wk
#	subset qtl object to region of interest
# LOD
lodA <- qtl$lod$A
lodA <- lodA[lodA$chr == 2,]
lodA <- lodA[lodA$pos > 101.7, ]
lodA <- lodA[lodA$pos < 113.7, ]
# Coef
coefA <- qtl$coef$A #matrix
coefA <- as.data.frame(coefA)
coefA <- coefA[rownames(lodA),]
# replace 
qtl$lod$A <- lodA
qtl$coef$A <- coefA
#	plot allele effect 
qtl$coef$A[abs(qtl$coef$A) > 2 ] = 0 
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/coef.Alb10.chr2.bayesian.int.pdf", width = 10.0, height = 7.5)
coefplot(qtl, chr = 2, main = "Chr 2 Alb10WK_log Allele Effect Plot @ bayesian interval")
dev.off()

## Ploting average founder effect of each strain
#	subset out coefficent data
data <- qtl$coef$A
colnames(data)[1] <- "A"
data$sex <- NULL
data$creat10wk <- NULL
#	Center the coefficent values
data[,2:ncol(data)] <- data[,2:ncol(data)] + data[,1]
data <- data - rowMeans(data)
#	reshpe melt dataframe for ggplot
Founder_names <- c("A", "B", "C","D","E","F","G","H")
ggdata <- melt(data, variable.name = "Founders", value.name = "Effect")
#	Get standard error
ggdata_SE <- summarySEwithin(ggdata, measurevar = "Effect", withinvars = "Founders")

ggplot(ggdata, aes(Founders, Effect)) +
		geom_point()+
		geom_errorbar(aes(ymin = Effecct - se, ymax = Effect + se), width = 0.5)


pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Avg_founder_effect_QTLAlb6wk_BayInt_Chr2.pdf", width = 10.0, height = 7.5)
ggplot(ggdata_SE, aes(Founders, Effect)) +
	geom_point()+
	geom_errorbar(aes(ymin = Effect - se, ymax = Effect + se), width = 0.5) +
	labs( title = "Average founder effect within bayesian interval of Alb10wk QTL Chr 2", x = "DO Founders", y = "Founder Effect") +
	theme( plot.title = element_text(hjust = 0.5))
dev.off()



