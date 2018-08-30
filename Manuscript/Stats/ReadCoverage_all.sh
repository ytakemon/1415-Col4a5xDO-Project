# Yuka Takemon
# 08/29/18
# Run each sample with ReadCoverage_all.sh

cd /projects/korstanje-lab/ytakemon/Col4a5xDO/Scripts
while read sample
do
    qsub -v sample=${sample} ReadCoverage_one.sh
done < /projects/korstanje-lab/ytakemon/Col4a5xDO/civet_run/sub_Name_list.txt
