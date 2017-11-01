# Install JunctionSeq
# Documentation: https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("JunctionSeq")

# Install QoRT
# Documentation: http://hartleys.github.io/QoRTs/index.html
install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",
                   repos=NULL,
                   type="source")
