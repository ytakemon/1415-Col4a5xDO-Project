#QC_GeneSeek README
This directory contains scripts that were used to QC GeneSeek genome probability data. These are my first ever scripts so documentations aren't very good.

#Description of Scripts:
###[Best.1415.genprob.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/QC_GeneSeek/Best.1415.genprob.R)
* Identified best GeneSeek files, Gillian's files, and RNA-seq reconstruction samples to use.
* Compiled into subsets of best samples 
* Established best.genoprobs.192 samples used throughout analysis.

###[Best.genoprob.192.kinship.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/QC_GeneSeek/Best.genoprob.192.kinship.R)
* Final kinship check of best.genoprobs.192 before using for analysis.

###[Gillian.compare.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/QC_GeneSeek/Gillian.compare.R)
* Comparing Gillian's samples with our RNA-seq reconstruction. 
* Data matched.
* Conclusion, our DO samples were swapped with her B6 samples.

###[reconst.v.geneseek.cor.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/QC_GeneSeek/reconst.v.geneseek.cor.R)
* Comparing all reconstruction data with GeneSeek genome probability data.
* Identified problem samples received from GeneSeek.
