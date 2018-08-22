cd /projects/korstanje-lab/ytakemon/Col4a5xDO/Scripts
while read sample
do
    qsub -v sample=${sample} Samtools_sort_one.sh
done < /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/Name_list.txt

# OOPS!
#while read sample
#do
#   rm /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/*${sample}/${sample}.bowtie.sorted.bam
#done < /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/Name_list.txt
