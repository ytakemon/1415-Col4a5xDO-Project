# Install QoRT Java file from: http://hartleys.github.io/QoRTs/QoRTs-STABLE.jar
#
# Documented usage:
# java -jar /path/to/jarfile/QoRTs.jar QC input.bam anno.gtf.gz /output/dir/
# Command line help:
# Additional options and syntax information for the main QC java utility can be found using the command:
java -jar /path/to/jarfile/QoRTs.jar QC --man
#Options and information about other sub-utilities within the java package can be found using the command:
java -jar /path/to/jarfile/QoRTs.jar --man
#And for each sub-utility:
java -jar /path/to/jarfile/QoRTs.jar utilname --man

# Bam files must be sorted. How to check: https://www.biostars.org/p/5256/
# % samtools view -H 5_110118_FC62VT6AAXX-hg18-unsort.bam
# @HD    VN:1.0    SO:unsorted
# % samtools view -H 5_110118_FC62VT6AAXX-hg18-sort.bam
# @HD    VN:1.0    SO:coordinate
