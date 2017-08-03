# Correlation README  
This directory contains correlation of TPMs calcualted from from RNA-Seq civet output.

# Description of Scripts:
### [C2.GFR.GE.cor.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/Correlation/C2.GFR.GE.cor.R)
* Started as correlation between GFR and all tpm.
* Restarted with PCA QC step in Expression correlation with phenotype.R
* Should be moved to scrap

### [Col12345comparison.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/Correlation/Col12345comparison.R)
* Correlation and ANOVA between Type IV Collagen expression (a1-4).
* Correlation and ANOVA between Type IV Collagen and phenotype.
* Plot by genotype and sex.

### [Expression correlation with phenotype.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/Correlation/Expression%20correlation%20with%20phenotype.R)
* PCA QC step to make sure there are no unknown factors skewing our data.
* Only identified sex, and a couple outliers that were due to abnormal expression of X and Y specific genes.
* Correlation and ANOVA between TPM and phenotypes ( GFR, ACR, ALB at all time points).
* Contains BioMart accesss to convert ENSEMBL mouse ID to more human readable MGI IDs.

### [Fmn1_gremlin_comparison.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/Correlation/Fmn1_gremlin_comparison.R)
* Correlation and ANOVA between Fmn1 and Gremlin.

### [Rfx3_Fmn1_correlation.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/Correlation/Rfx3_Fmn1_correlation.R)
* Correlation and ANOVA between Rfx3 and GFR.
* Correlation and ANOVA between FMN1 and ACR at 10wks.

###
[Tgm2_correlation.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/RNA_seq/Correlation/Tgm2_correlation.R)
* Correlation and ANOVA between Tgm2 and both GFR and ACR.
