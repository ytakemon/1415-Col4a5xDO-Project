#PBS -l nodes=1:ppn=3,walltime=12:00:00

# Yuka Takemon
# 08/29/18
# Sort and calcuate depth

module load samtools/1.8
cd /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/*${sample}

# sort bam file
samtools sort ${sample}.bowtie.bam > ${sample}.bowtie.sorted.bam

# get depth using -a to also gather 0 reads
samtools depth -a ${sample}.bowtie.sorted.bam > ${sample}.bowtie.sorted.bam.samtools.depth

# Get rid of sorted bam, takes up too much space of cluster
rm ${sample}.bowtie.sorted.bam
