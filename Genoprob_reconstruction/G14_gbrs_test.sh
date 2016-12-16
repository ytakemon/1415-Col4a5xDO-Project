#reconstruct genoprob using Civet pipeline EMASE data

#!/bin/bash -l
#PBS -l nodes=1:ppn=4,walltime=24:00:00

#load modules
module load Anaconda/2.4.0
source activate gbrs
export GBRS_DATA=/home/ytakemon/GBRS_DATA

#set directories
civet_emase_dir=/hpcdata/ytakemon/Col4a5xDO/civet_run
GBRS_DATA=/home/ytakemon/GBRS_DATA/R84-REL1505

gbrs reconstruct -e ${civet_emase_dir}/*1415-0674*/gbrs.quantified.multiway.genes.tpm \
		 -t ${GBRS_DATA}/tranprob.DO.G14.F.npz \
		 -x ${GBRS_DATA}/avecs.npz \
 		 -g ${GBRS_DATA}/ref.gene_ids.ordered.npz  

