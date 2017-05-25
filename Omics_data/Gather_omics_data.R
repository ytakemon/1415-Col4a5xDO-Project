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

# A: 4.5 weeks
# B: 6.0 weeks
# C: 7.5 weeks
library(stringr)
setwd("/Users/ytakemon/GitHub/1415-Col4a5xDO-Project")

Sieve20 <- read.delim("Omics_data/Data/Muckova/AS0_C0_Sieve2.0.txt", sep = "\t", header = TRUE)
Sieve20 <- Sieve20[,c("Protein.ID", "Description", "Peptides", "Frames",
                      "Hits", "Ratio.AS.0.C.0", "StDev", "pValue","Fraction.1D_2D")]

Sieve20$Gene_name <- str_split_fixed(Sieve20$Description, "GN=", 2)
