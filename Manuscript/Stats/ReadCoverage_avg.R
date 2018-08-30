# Yuka Takemon
# 08/29/18
# Get average of read coverage after ReadCoverage_one is done
wd <- "/projects/korstanje-lab/ytakemon/Col4a5xDO/"
setwd(wd)

# read data
df <-  read.delim("./civet_run/AvgReadCoverage.txt", sep ="\t", stringsAsFactors = F, header = F)
mean(df[,2])
#> mean(df[,2])
#[1] 37.663
