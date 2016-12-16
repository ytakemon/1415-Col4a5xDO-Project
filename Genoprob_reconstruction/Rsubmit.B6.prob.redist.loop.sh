#this script can be run directly on cadillac via bash scriptname.sh

#!/bin/bash
cd /hpcdata/ytakemon/Col4a5xDO/Scripts

while read I
do
qsub -v I=${I} B6.prob.redistribution.sh

done < /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/DOQTL.GM.interpol.reconstruction/X1415.names.only.txt
