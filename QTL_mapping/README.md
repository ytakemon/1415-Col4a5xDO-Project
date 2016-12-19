#QTL mapping README
This directory contains all scripts used for QTL analysis for Albuminuria and GFR.

#Description of Scripts:
### Ha_QTL_ACR_Alb.R: 
* Using best.genoprobs.192 data.
* Contains Ha QTL analysis for ACR and Albumin with sex (or sex and creatinine) as additive covariates. Additionally permutations, thresholding and plotting are included. 
* Added script to make QTL plot excluding the X-chromosome for publication.
* Added script to make allele effect plot to focus on region of interest rather than plotting whole chromosome.

### Ha_QTL_compare_GFR_calc.R
* Using best.genoprobs.192 data.
* Contains Ha QTL analysis for GFR models C2, 2xC1, C1, and PL with sex as an additive covariate.
* Contains premutations, thresholding, and plotting.
* Added script to plot QTL without X-chromosome.

### Ha_QTL_GFR_C2.R
* Determined best.genoprobs.192 data in this script.
* Contains Ha QTL analysis for C2 GFR comparing different subset of samples to determine most optimal genoprob sample size given GeneSeek mistake.
* Subsets include 192 samples (157 GeneSeek + 24 RNA-seq reconstruct + 11 Gillian samples), only 157 GeneSeek samples (confirmed to be accurate), and 182 GeneSeek Samples ( includes GeneSeek mistakes).

### Ha_Sex_specific_QTL_GFR_Alb.perms.R
* Uses best.genoprobs.data and Ha QTL data previously calculated for GFR and Alb.
* Samples wer subsetted to females and males for single sex only analysis.
* INcludes permutaiton, thresholding, and plotting.

### Hf_Hi_QTL_sex_intcovar_Alb.R
* Using best.genoprobs.192 data
* Contains Hf QTL analysis of Albumin to determine sex effect using sex as an interactive covariate, while keeping creatinine as an additive covariate.
* Grabs permutation results created by Hf_QTL_sex_intcovar_Alb_perm.R
* Creates threshold and plots QTL.
* Added Hi QTL analysis using Hf QTL object as a shell, only lod values are accurate to Hi analysis.
* Contains Hi QTL permutation created by Hi_QTL_Alb_perm.R, thresholding, and qtl plotting.

### Hf_Hi_QTL_sex_intcovar_GFR.R
* Using best.genoprobs.192 data.
* Contains Hf QTL analysis of C2 GFR to determine sex effect using sex as an interactive covariate.
* Grabs permutation results created by Hf_QTL_sex_intcovar_GFR_perm.R
* Creates threshold and plots QTL.
* Added Hi QTL analysis using Hf QTL object as a shell, only lod values are accurate to Hi analysis.
* Contains Hi QTL permutation created by Hi_QTL_GFR_perm.R, thresholding, and qtl plotting.

### Hf_QTL_sex_intcovar_Alb_perm.R
* Manual loop for Albumin Hf QTL permutation.
* Use in conjunction with Hf_QTL_sex_intcovar_Alb.R

### Hf_QTL_sex_intcovar_GFR_perm.R
* Manual loop for GFR Hf QTL permutation 
* Use in conjunction with Hf_QTL_sex_intcovar_GFR.R

### Hi_QTL_Alb_perm.R
* Manual loop for Albumin Hi QTL permuation.
* Use Hf QTL as shell only lod values are accurate to Hi.
* Use in conjunction with Hi_QTL_Alb.R.

### Hi_QTL_GFR_perm.R
* Manual loop for GFR Hi QTL permutaiton.
* Use Hf QTL as shell only lod values are accurate to Hi.
* Use in conjunction with Hi_QTL_GFR.R.
