
dir=/data
# Set the reference fasta
#phase2 reference
ref=hs37d5.fa

# Create Variable for input file
INPUT1=NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam
#choose how to log time
Time=/usr/bin/time

./freebayes/bin/freebayes \
-f ${dir}/${ref} \
-v ${dir}/${INPUT1}.freebayes.vcf \
${dir}/${INPUT1} \
> freebayesADAM.out 2> freebayes.err

#freebayes reference output input