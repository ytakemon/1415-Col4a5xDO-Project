library(stringr)
setwd("/projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/")
list <- list.files(pattern = "fastq.gz")
list <- str_sub(list,15,)

write.table(list, "./Name_list.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
