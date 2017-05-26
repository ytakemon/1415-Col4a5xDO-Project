# Gather public kidney datasets to check if we have genes that overlap.
# The data sets collected are the following:
##  HKUPP and UPDB refered in the Papadopoulos et al. paper.
##  Proteomics dataset from Col4a3 KO mice in Muckova paper.
##  Rinchen proteomics datasets: Isolated rat podocytes from induced kidney damage, human
##  podocyte cell line with induced damage.

# Extensive records should be kept to make sure we interpret data correctly.

# Preclinical Alterations in the Serum of COL(IV)A3–/– Mice as Early Biomarkers of Alport Syndrome
# Petra Muckova†‡, Sindy Wendler†, Diana Rubel§, Rita Büchler†, Mandy Alert†, Oliver Gross§, and Heidrun Rhode*†
# † Institute of Biochemistry I, University Hospital Jena, Nonnenplan 2-4, 07740 Jena, Germany
# ‡ Clinic of Neurology, University Hospital Jena, Erlanger Allee 101, 07740 Jena, Germany
# § Department of Nephrology and Rheumatology, University Medicine Göttingen, Robert-Koch Str. 40, 37075 Göttingen, Germany

# Understanding Sieve analysis by thermo fisher (http://tools.thermofisher.com/content/sfs/manuals/Man-XCALI-97186-SIEVE-110-User-ManXCALI97186-B-EN.pdf)
## SIEVE is an automated software package for the label-free, semi-quantitative differential
## expression of proteins and peptides. It performs comparative analyses of sample populations
## by comparing the raw spectral information from LC/MS analyses of control (healthy) and
## disease (treatment) samples to determine if there are changes among the two sample sets
## which may indicate differential protein expression.

# Values in from Sieve 2.0 data provides differential expression of proteins expressed by the ratio of
# AS over the control. I would need to talk to Dorothy at mass spec to make sure I am looking at the data correctly

# A: 4.5 weeks
# B: 6.0 weeks
# C: 7.5 weeks
library(stringr)
setwd("/Users/ytakemon/GitHub/1415-Col4a5xDO-Project")

Sieve20 <- read.delim("Omics_data/Data/Muckova/AS0_C0_Sieve2.0.txt", sep = "\t", header = TRUE)
Sieve20 <- Sieve20[,c("Protein.ID", "Description", "Peptides", "Frames",
                      "Hits", "Ratio.AS.0.C.0", "StDev", "pValue","Fraction.1D_2D")]

# Extract gene name from Description
Sieve20$Gene_name <- str_split_fixed(Sieve20$Description, "GN=", 2)[,2]
Sieve20$Gene_name <-str_split_fixed(Sieve20$Gene_name, " PE=", 2)[,1]

# Clean up protein ID (eg: SP://A2A891-5, A2A891 is the Uniprot protein ID, A2A891-5 identified isoform 5 of protein)
Sieve20$UniprotID <- str_split_fixed(Sieve20$Protein.ID, "SP://", 2)[,2]

#Checking how many proetins are significant:
# > dim(Sieve20)
# [1] 7402   11
# > dim(Sieve20[Sieve20$pValue < 0.05,])
# [1] 7342   11
# > dim(Sieve20[Sieve20$pValue < 0.001,])
# [1] 5466   11
# That's about a fifth of all proteins ... est. 25,000 proteins in mice

#Checking out Sieve2.1 data
Sieve21 <- read.delim("Omics_data/Data/Muckova/AS0_C0_Sieve2.1.txt", sep = "\t", header = TRUE)
Sieve21 <- Sieve21[,c("Protein.ID", "Description", "Peptides", "Frames",
                      "Hits", "Ratio.AS.0.C.0", "StDev", "pValue","Fraction.1D_2D")]

# Extract gene name from Description
Sieve21$Gene_name <- str_split_fixed(Sieve21$Description, "GN=", 2)[,2]
Sieve21$Gene_name <-str_split_fixed(Sieve21$Gene_name, " PE=", 2)[,1]

# Clean up protein ID (eg: SP://A2A891-5, A2A891 is the Uniprot protein ID, A2A891-5 identified isoform 5 of protein)
Sieve21$UniprotID <- str_split_fixed(Sieve21$Protein.ID, "SP://", 2)[,2]

# > dim(Sieve21)
# [1] 4103   11
# > dim(Sieve21[Sieve21$pValue < 0.05,])
# [1] 2603   11
# > dim(Sieve21[Sieve21$pValue < 0.01,])
# [1] 1547   11
# > dim(Sieve21[Sieve21$pValue < 0.001,])
# [1] 870  11
# Much more reasonable... First look at the data suggests some corrections were made, will have to
# go through literature to figure out what they did or if its a software update.
# Since I dont have the data to run a permutation test, I will use 3 thresholds 0.05, 0.01, and 0.001.

Sieve21_05 <- Sieve21[Sieve21$pValue < 0.05,]
Sieve21_01 <- Sieve21[Sieve21$pValue < 0.01,]
Sieve21_001 <- Sieve21[Sieve21$pValue < 0.001,]
