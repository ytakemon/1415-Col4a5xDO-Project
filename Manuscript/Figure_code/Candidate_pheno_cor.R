# Yuka Takemon
# 08/29/18
# Create correlation between phenotype and candiate genes
library(tidyverse)

# set directory
wd <- "/Users/ytakemon/Dropbox/Col4a5/Data"
setwd(wd)

#	load files
load("./Col4a5xDO_192data_YT.Rdata")

# function to get genes
# Get_gene()
# A function to extract data for one gene into a tibble, borrowed from Gary
Get_gene <- function(GeneSymbol, phenotype){

  # Get animal information
  MouseID <- Pheno$MouseID
  Sex <- factor(as.character(Pheno$Sex), levels=c("F","M"), labels=c("F","M"))


  # get gene ID from symbol
  ensID <- AnnotGenes[AnnotGenes$external_gene_name == GeneSymbol,"ensembl_gene_id"]

  # build data frame
  outputDF <- data_frame(MouseID = MouseID,
                         Sex = Sex)

  # Add expression data
  if(length(ensID) > 0){
    outputDF <- mutate(outputDF,
      RNA = mRNAexpr[,ensID])
  }

  # Add phenotype
  if(!missing(phenotype)){

    if( phenotype == "ACR"){
      outputDF <- bind_cols(outputDF, Pheno[,c("ACR6WK","ACR10WK","ACR15WK")])
    } else if( phenotype == "GFR"){
      outputDF <- mutate(outputDF,
        C2_log = Pheno$C2_log)
    }
  }

  # return data
  return(outputDF)
}
# #test
# Get_gene("Gnai3", "ACR")
# Get_gene("Cow")      # a non-existing gene

# Correlation between Rfx3 and GFR phenotypes
Rfx <- Get_gene("Rfx3", "GFR")

# Need to get allele info at loci for each sample
GetAllele <- function(df, gene){
  # Get gene annotation
  position <- AnnotGenes[AnnotGenes$external_gene_name == gene,]

  # if it doens't exist
  if(nrow(position) == 0 ){
    print(paste("Gene:", gene, "not found"))
    return(df)
  }

  sub_snp <- Snps %>% filter((chr == position$chromosome_name) &
                             (pos > position$start_position * 1e-6) &
                             (pos > position$end_position * 1e-6))

  # Get marker information and figure out most frequent allele for each sample
  # Gather markers in genomic region of region
  sub_probs <- Genoprobs[,,colnames(Genoprobs[,1,]) %in% sub_snp$marker]
  # Which allele is most frequent for each marker in this region?
  x <- NULL
  for(i in 1:dim(sub_probs)[3]){
    temp <- sub_probs[,,i]
    y <- NULL
    for(j in 1:nrow(sub_probs)){
      allele <- names(which(temp[j,] > 0.4))
      if(length(allele) > 1){
        allele <- paste0(allele[1], allele[2])
      }
      y <- c(y, allele)
    }
    x <- cbind(x,y)
  }
  x <- t(as.data.frame(x, stringsAsFactors = FALSE))
  colnames(x) <- rownames(sub_probs)
  rownames(x) <- colnames(sub_probs[1,,])
  # Which allele is most dominant in this region for each sample?
  Dom_Allele <- NULL
  for(i in 1:ncol(x)){
    dom <- names(which(table(x[,i]) == max(table(x[,i]))))
    if(length(dom) > 1){
      dom <- paste0(dom[1], dom[2])
    }
    Dom_Allele <- c(Dom_Allele, dom)
  }
  df$Allele <- Dom_Allele
  return(df)
}

# test
#GetAllele(Rfx, "Rfx3")
#GetAllele(Rfx, "Cow")

#
Rfx <- GetAllele(Rfx, "Rfx3")
Rfx$FounderEffect <- NA
Rfx$FounderEffect[which(Rfx$Allele %in% c("A","D","C"))] <- "Low"
Rfx$FounderEffect[which(Rfx$Allele %in% c("G","B"))] <- "High"

Rfx %>% filter(Sex == "M") %>%
ggplot(., aes(x = RNA, y = C2_log, colour = FounderEffect))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, colour = "black")+
  xlab("Rank transformed Rfx3 TPM")+
  ylab("Log-transformed GFR")
