# Gather public kidney datasets to check if we have genes that overlap.
# The data sets collected are the following:
##  HKUPP and UPDB refered in the Papadopoulos et al. paper.
##  Proteomics dataset from Col4a3 KO mice in Muckova paper.
##  Rinchen proteomics datasets: Isolated rat podocytes from induced kidney damage, human
##  podocyte cell line with induced damage.

# Extensive records should be kept to make sure we interpret data correctly.

# YAP-mediated mechanotransduction determines the podocyte’s response to damage
# Markus M. Rinschen,1,2,3,4 Florian Grahammer,5,6* Ann-Kathrin Hoppe,1,2*
# Priyanka Kohli,1,2,3 Henning Hagmann,1,2 Oliver Kretz,5,6,7 Sabine Bertsch,1,2
# Martin Höhne,1,2,3,4 Heike Göbel,8 Malte P. Bartram,1,2 Rajesh Kumar Gandhirajan,1
# Marcus Krüger,2,3 Paul-Thomas Brinkkoetter,1,2 Tobias B. Huber,5,6,9,10
# Martin Kann,1,2 Sara A. Wickström,3,11
# Thomas Benzing,1,2,3,4† Bernhard Schermer1,2,3,4†

# Data information provided by author:
# This data demonstrates the ttest differences of logarithmized (log2) expression
# values of every identified protein and the respective pvalue.
# This matrix is based on the proteingroups.txt output of MaxQuant software.
# Each column represents a different animal.
# For details and more information, please check the raw data deposited at PRIDE
# and the methods section of the manuscript.

# Data retreived from supplement:
# aaf8165_Data file S1.xlsx <- Proteomic analysis of PAN effect in rat glomeruli in Excel format.
# aaf8165_Data file S2.txt <- Proteomic analysis of PAN effect in rat glomeruli in text format.
# aaf8165_Data file S3.xlsx <- Proteomic analysis of PAN effect in human podocytes in Excel format.
# aaf8165_Data file S4.txt <- Proteomic analysis of PAN effect in human podocytes in text format.

# Notes:
# Dont use excel sheet, some gene names were converted to dates

#Rat glomeruli
library(stringr)
setwd("/Users/ytakemon/GitHub/1415-Col4a5xDO-Project")

Rat_glom <- read.delim("Omics_data/Data/Rinschen/aaf8165_Data file S2.txt", sep= "\t")
Rat_glom <- Rat_glom[-c(1:3),]
names(Rat_glom)[1:12] <- c("con1", "con2", "con3", "con4",
                            "2d1", "2d2", "2d3", "2d4",
                            "4d1", "4d2", "4d3", "4d4")
Rat_glom$Gene.names <- as.character(Rat_glom$Gene.names)
Rat_glom$Gene.names[Rat_glom$Gene.names == ""] <- NA
Rat_glom <- Rat_glom[complete.cases(Rat_glom$Gene.names),]
Rat_glom$Gene_name <- str_split_fixed(Rat_glom$Gene.names, ";", 2)[,1]
Rat_glom$ProteinID <- str_split_fixed(Rat_glom$Protein.IDs, ";", 2)[,1]

#Need to redo t-test between the 3 groups because file does not give pval, only test stat (ugh... *eye roll*)
