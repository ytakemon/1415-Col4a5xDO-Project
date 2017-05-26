# I need to collect all the gene names, alias, and convert them all to Ensembl IDs
# for stream lined query.
# https://www.biostars.org/p/14971/#212906
# Also check out WormHole downloads and BioMart

# Download mouse annotation from bioconductior
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
library(RSQLite)

queryGeneNames <- c("Krt6b")
dbCon <- org.Mm.eg_dbconn()
sqlQuery <- "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id"
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
result <- aliasSymbol[which(aliasSymbol[,2] %in% queryGeneNames),5]

#Should make this into an independent script to run in ther Rscripts. 
