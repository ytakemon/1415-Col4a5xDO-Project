#PBS -l nodes=1:ppn=3,walltime=12:00:00
module load samtools/1.8
cd /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/*${sample}

# sort bam file
samtools sort ${sample}.bowtie.bam > ${sample}.bowtie.sorted.bam

# get depth using -a to also gather 0 reads
samtools depth -a ${sample}.bowtie.sorted.bam > ${sample}.bowtie.sorted.bam.samtools.depth
