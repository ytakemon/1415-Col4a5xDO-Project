# eQTL README
This directory contains scripts for eQTL analysis.

# Description of scripts:
### [Allele_effect_trans_eQTL.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/eQTL/Allele_effect_trans_eQTL.R)
* Allele effect plot for eQTL identified to have significant trans-eQTL at Fmn1.

### [eQTL.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/eQTL/eQTL.R)
* Contains multiple single eQTL analysis of initial genes of interest

### [eQTL_all.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/eQTL/eQTL_all.R)
* Contains a loop script to analyze all RNA-seq TPMs for eQTL.
* Creates and saves qtl (scanone) object for every gene expressed,
* Creates and saves pdf of qtl manhattan plots. 

### [trans_eQTL.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/eQTL/trans_eQTL.R)
* contains a loop script to identify LOD score at FMN1 (region of interst).
* Then subsetted into genes that surpass 3 levels of LOD significance. 
