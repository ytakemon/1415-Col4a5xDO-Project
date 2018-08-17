#!/bin/bash -l
#PBS -l nodes=1:ppn=3,walltime=12:00:00
cd /hpcdata/ytakemon/Col4a5xDO/Scripts/
module load R/3.4.4

Rscript ReadCoverage.R ${sample}
