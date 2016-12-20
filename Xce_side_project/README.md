#Xce_side_project
This directory contains scripts that pertains to our X controlling element (Xce) side project. We know that various Xce allele 
determines the likelihood in which X chromosome becomes silenced in females. We want varify Xce allele in our Col4a5xDO F1 animals
and analyse its effect on severity of renal phenotype. 

# Description of scripts:
### [QTL_Alb_Xce.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Xce_side_project/QTL_Alb_Xce.R)
* Trial using Xce as another additive covariate.

### [Xce.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Xce_side_project/Xce.R)
* Determining Xce allele of each animal
* Determining allele of Xist expression. 
* Comparing Xce to Xist expression allele to make sure that Xist is expressed from the silenced gentype. 
* Correlation and ANOVA test for Xce and phenotype.

### [sex.determination.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Xce_side_project/sex.determination.R)
* Determining sex of each sample by expression of X and Y specific markers to address issues identified in PCA analysis in RNA_seq/Correlation/Expression correlation with phenotype.R 
