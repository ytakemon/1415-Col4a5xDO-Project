##Xce and Xist allele for Col4a5xDO
#	Yuka Takemon
#	Created: 11/09/16
#	Modified last: 21/10/20

##Intro:
#	We know that X chromosome silencing is driven by a region called the Xce (X controlling element), which is a cis-regulatory gene.
#	Xist, 5' to Xce, is expressed by the inactive X-chromosome and facilitates the spread of inactivation of that chromosome.
#	In Mus musculus there are 4 known Xce alleles: Xce^a, Xce^b, Xce^c, and Xce^e, which effect the probability that the X-chromosome will be
#	expressed. The liklihood is as follows: Xce^a < Xce^e < Xce^b < Xce^c
#	We would like to perform ANOVA test to see if the allele of Xce allele has an effect on the phenotype on females as they are Hets for the Col4a5
#	mutation located on the B6 X-chromosome.

##Xce allels of 8 founders:
#	Xce^a : AJ, 129
#	Xce^b : WSB, B6, NOD, NZO
#	Xce^c : CAST
#	Xce^e :	PWK

##Expectation of phenotype severity:
#	Xce^a x Xce^b: will have more severe phenotype as the mutation is on the Xce^b chromosome.
#	Xce^b x Xce^b: will have equal distribution of phenotype as the probability is equal.
#	Xce^c x Xce^b: will have less severe phenotype as the dominant allele does not carry the mutaiton.
#	Xce^e x Xce^b: will have more severe phenotype as the Xce^e is less dominanat that the allele carrying the mutation.

##Xce is not yet on Ensemble but there are a few papers that have located a specific locus that suggests to the the Xce alelle. (Calaway et al., 2013)
#	Using the two markers: JAX00716733 and JAX00182728 we can distinguish between the four alleles.
#	Xce allele JAX00716733 JAX00182728
#	Xce^a 		C 			C
#	Xce^b 		C 			T
#	Xce^c 		G 			G
#	Xce^e 		G 			T

##Xist vs Xce
#	We can check the effects of Xce by identifying the gene expression of Xist gene and which gene it was derived from as the EMASE is a able to account
#	for this. Xist should be expressed form the X chromosoem that is being suppressed (the weaker Xce allele)

##Workflow:
#	First identify the snp region and assign specific alleles at the two markers using Genoprobs data (I think)
#	Then using the EMSE output identify founders that contirubted to the Xist expression
#	Compare Xce and Xist
#	Compare Xce to phenotypes (GFR & ACR)
#	Compare Xist to phenotypes (GFR & ACR)

##Xce interval
#	x:102826982-102990039

#setup:
library(pacman)
p_load(DOQTL, tidyverse, ggsci)

dir <- "/projects/marralab/ytakemon_prj/Col4a5/"
#load sample
load(paste0(dir, "Data/consolidated/best.genoprobs.192.Rdata"))
load(paste0(dir, "Data/consolidated/GM_snps.Rdata"))
load(paste0(dir,"Data/consolidated/Phenotype/Xce_allele.Rdata"))
pheno <- read.delim(paste0(dir,"Data/consolidated/Phenotype/1415_master_pheno.txt"), sep = "\t", header = TRUE)

# Clean up pheno data ------------------------------------------------
rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names

# Add Xce allele info
pheno <- pheno[rownames(Xce_allele),] #subset pheno to match only females
pheno$Winner_origin <- Xce_allele$Winner_origin
pheno$Xce_DO <- Xce_allele$Xce_DO
pheno$Xce_winner <- Xce_allele$Xce_winner

# Get ACR values
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$C2_log <- log(pheno$C2)
pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR10WK_log <- log(pheno$ACR10WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs

# Figure 3 plot
ggdata <- pheno %>%
  mutate(sample = rownames(pheno)) %>%
  select(sample, Sex, Xce_DO, ACR6WK_log, ACR10WK_log, ACR15WK_log, C2_log) %>%
  rename(log_C2 = C2_log) %>%
  pivot_longer(-c(sample,Sex,Xce_DO, log_C2), names_to = "ACR_month", values_to = "log_ACR_value") %>%
  mutate(Age_in_Months = case_when(
            .$ACR_month == "ACR6WK_log" ~ "6",
            .$ACR_month == "ACR10WK_log" ~ "10",
            .$ACR_month == "ACR15WK_log" ~ "15",
            TRUE ~ "NA")) %>%
  mutate(Age_in_Months = fct_relevel(Age_in_Months, "6","10","15")) %>%
  select(sample, Sex, Age_in_Months, Xce_DO, log_C2, log_ACR_value)

# Plot figure 3
pdf(paste0(dir,"Results/Xce_allele_ACR.pdf"), height = 4, width = 6)
ggplot(ggdata, aes(x = Age_in_Months, y = log_ACR_value, fill = Xce_DO))+
  geom_boxplot() +
  scale_fill_aaas() +
  labs(x = "Weeks of age",
       y = "log(ACR)",
       fill = "Xce allele") +
  theme_bw()+
  ylim(1,10)
dev.off()

# Get sample numbers for Figure 3
ggdata %>% drop_na() %>% count(Age_in_Months)
ggdata %>% drop_na() %>% count(Age_in_Months, Xce_DO)

# export data for sharing:
write_csv(ggdata, paste0(dir,"Results/Sample_Xce_info_with_phenotype.csv"))

# get mean and SD for GFR
ggdata %>% select(sample, Sex, Xce_DO, log_C2) %>% group_by(Sex, Xce_DO) %>%
  na.omit() %>%
  summarise(log_C2_mean = mean(log_C2),
            log_C2_sd = sd(log_C2))


# > fit <- aov(log_C2 ~ Xce_DO, df)
# > summary(fit)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Xce_DO        3  1.764  0.5880   3.553 0.0155 *
# Residuals   191 31.609  0.1655
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > TukeyHSD(fit)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = log_C2 ~ Xce_DO, data = df)
#
# $Xce_DO
#                   diff         lwr        upr     p adj
# Xce_b-Xce_a -0.0762419 -0.25020279 0.09771899 0.6678057
# Xce_c-Xce_a  0.1524694 -0.05568201 0.36062078 0.2321434
# Xce_e-Xce_a -0.2239204 -0.59586442 0.14802364 0.4038323
# Xce_c-Xce_b  0.2287113  0.01908161 0.43834097 0.0264456
# Xce_e-Xce_b -0.1476785 -0.52045182 0.22509486 0.7339746
# Xce_e-Xce_c -0.3763898 -0.76629306 0.01351351 0.0627897

pdf(paste0(dir,"Results/Xce_allele_GFR.pdf"), width = 6, height = 4)
ggplot(ggdata, aes(x = Xce_DO, y = log_C2, fill = Xce_DO))+
  geom_boxplot() +
  scale_fill_aaas() +
  labs(x = "Xce allele",
       y = "log(GFR)",
       fill = "Xce allele") +
  theme_bw()
dev.off()
