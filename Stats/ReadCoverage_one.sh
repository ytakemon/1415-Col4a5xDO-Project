#!/bin/bash -l
#PBS -l nodes=1:ppn=3,walltime=15:00:00

# Yuka Takemon
# 08/29/18
# Run ReadCoverage.R for each sample

cd /projects/korstanje-lab/ytakemon/Col4a5xDO/Scripts
module load R/3.4.4

Rscript ReadCoverage.R ${sample}
