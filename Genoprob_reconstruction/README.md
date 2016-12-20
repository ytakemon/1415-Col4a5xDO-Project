# Genoprob_reconstruciton README  
This directory contains scripts that were used to reconstruct genomes/genome probabilities using RNA-Seq from multiple parent (Col4a5xDO) mice.
References to the GBRS suit used were written by K.B.Choi and can be found [here](https://media.readthedocs.org/pdf/gbrs/latest/gbrs.pdf) and [documentation here](http://gbrs.readthedocs.io/en/latest/usage.html#to-use-gbrs-in-command-line).

# Description of Scripts:
###[submit_civetSE.sh](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/submit_civetSE.sh)
* Bash script.
* civet run for single end RNA-seq genome reconstruction.

###[G14_gbrs_test.sh](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/G14_gbrs_test.sh)
* Bash script.
* Genome reconstruction from RNA-seq data.
* Test run for one sample (1415-0674).
* Using civet pipleline EMASE output data.

###[PySubmit_subset_reconstruct.sh](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/PySubmit_subset_reconstruct.sh)
* Bash script.
* Genome reconstruction from RNA-seq data.
* Ran subset of samples that were identified to replace problematic GeneSeek data.
* Using civet pipeline EMASE output data.
* Written in for loop format.

###[Rsubmit.B6.prob.redist.loop.sh](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/Rsubmit.B6.prob.redist.loop.sh)
* Bash script.
* Loop script to run B6.prob.redistribution.sh.

###[B6.prob.redistribution.sh](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/B6.prob.redistribution.sh)
* Bash script.
* Calls resources on cadillac and runs Rscript B6.prob.redistribution.R

###[B6.prob.redistribution.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/B6.prob.redistribution.R)
* Rscript.
* Takes interpolated data from interpolate.64k2GM.R
* Redistributes 50% of the B6 probability contributed by the Col4a5 mutant (on B6 background) to focus probabily on the 8 founders.
* Divides 50% from B6 and divides everything by the sum of the row. 

###[interpolate.64k2GM.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/interpolate.64k2GM.R)
* Rscript.
* interpolates 64K grid made by genome reconstruction into GM_snp grid used by GeneSeek.

###[export-genoprob-tsv.sh](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/export-genoprob-tsv.sh)
* Bash script
* Takes genome reconstruct and interpolates up to standard 64K grid.

###[genoprobs_reconstruction.R](https://github.com/TheJacksonLaboratory/1415-Col4a5xDO-Project/blob/master/Genoprob_reconstruction/genoprobs_reconstruction.R)
* Rscript.
* Combined reconstruct with GeneSeek genoprob.
* Includes kinship QC.
