#!/bin/bash -l
#PBS -l  nodes=1:ppn=4,walltime=02:00:00 
cd /hpcdata/ytakemon/Col4a5xDO/Scripts/Rscripts

module load R/3.3.1

R --no-save --args ${I} < B6.prob.redistribution.R > B6.prob.redistribution.${I}.Rout
