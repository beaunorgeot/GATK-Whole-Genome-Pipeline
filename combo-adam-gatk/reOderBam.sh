cd /data/NAparquet/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bqsr.bam
#cat the parts into a single bam mv the file to the main dir
samtools cat -o combinedNA12878.bam part-r-*
mv combinedNA12878.bam /data/


#Perform reorder in /data (or the dir that has the bad bam and the new reference)
cd /data

java -jar ~/picard.jar ReorderSam \
    I=combinedNA12878.bam \
    O=kayCombinedNA12878.bam \
    REFERENCE=hs37d5.fa