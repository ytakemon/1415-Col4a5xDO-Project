#Fit to 64k no Y no MT grid

#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=1:00:00

module load Anaconda/2.4.0
source activate gbrs

sample_list=/hpcdata/ytakemon/Col4a5xDO/Scripts/civet_list.txt

while read sample ;
do 
export-genoprob-file    -v \
                       	-i /hpcdata/ytakemon/Col4a5xDO/civet_run/${sample}/gbrs.interpolated.genoprobs.npz \
                        -s A,B,C,D,E,F,G,H \
                        -g /home/ytakemon/GBRS_DATA/R84-REL1505/ref.genome_grid.64k.noYnoMT.txt
done < $sample_list

