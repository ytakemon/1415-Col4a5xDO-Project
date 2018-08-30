# Yuka Takemon
# 08/29/18
# List all samples
library(stringr)
setwd("/projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/")
list <- list.files(pattern = "fastq.gz")
list <- str_sub(list,15,)

write.table(list, "./Name_list.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

# Get subset of samples that were not successfully run
library(stringr)
setwd("/projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/")
list <- read.delim("/projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/Name_list.txt",
  sep = "\t", header = F, stringsAsFactors = FALSE)
now <- read.delim("/projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/AvgReadCoverage.txt",
  sep = "\t", header = F, stringsAsFactors = FALSE)

list$name <- str_sub(list[,1], ,-17)
sublist <- list[!(list$name %in% now[,1]),]
sublist <- sublist[, -2]
write.table(sublist, "./sub_Name_list.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)




#set.seed(9269)
#ss <- sample(1:3,size=nrow(mtcars),replace=TRUE,prob=c(0.6,0.2,0.2))
#mycars <- setNames(split(mtcars,ss), c("train","test","cvr"))
