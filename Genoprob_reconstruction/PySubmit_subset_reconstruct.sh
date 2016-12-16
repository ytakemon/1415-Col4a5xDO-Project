#subset determined for genoprob reconstruction

#!/bin/bash

module load Anaconda/2.4.0
source activate gbrs

#Orignally 37 samples for reconstruction, but only 34 were sent to RNA-seq
sample_list = ( 1415-0170 1415-0176 1415-0198 1415-0201 1415-0212 1415-0226 1415-0235 1415-0248 1415-0258 1415-0259 1415-0261 1415-0266 1415-0302 1415-0318 1415-0324 1415-0329 1415-0340 1415-0347 1415-0348 1415-0351 1415-0372 1415-0429 1415-1014 1415-1029 1415-1035 1415-1043 1415-1051 1415-1054 1415-0256 1415-0194 1415-0674 1415-0457 1415-0807 )

#set directory
civet_output_dir = /hpcdata/ytakemon/Col4a5xDO/civet_run

for sample in ${sample_list}
gbrs reconstruct -e ${civet_outputdir}/*${sample}*/gbrs.quantified.multiway.genes.tpm \
		 -t ${civet_outputdir}/*${sample}*/



