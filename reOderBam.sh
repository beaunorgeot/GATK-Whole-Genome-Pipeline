
#Perform this in /data (or the dir that has the bad bam and the reference)

java -jar ~/picard.jar ReorderSam \
    I=combinedNA12878.bam \
    O=kayCombinedNA12878.bam \
    REFERENCE=hs37d5.fa