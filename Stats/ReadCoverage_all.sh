cd /projects/korstanje-lab/ytakemon/Col4a5xDO/Scripts
while read sample
do
    qsub -v sample=${sample} ReadCoverage_one.sh
done < /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/Name_list.txt
