#! /bin/bash
# This scipt is meant to take a file that has been preprocessed using ADAM and run # it through GATK's Haplotype variant caller

# The input is assumed to have been
# concatened using samtools cat,
# karyotypically ordered using my reOrderBam.sh script which takes ~12hrs on 234GB file


# NOTE: I SHOULD COMBINE SAMTOOLS CAT, REORDERBAM AND THIS INTO 1 SCRIPT "TRANSITION SCRIPT" and pass that to HAP-script
# The input when I compete this would then be combined.bam, not kayCombined.bam

set -e
set -x
set -o pipefail

cd ~

# Set $RAM to use for all java processes
RAM=-Xmx200g
# Set number of threads to the number of cores/machine for speed optimization
THREADS=32
# Set the dir for the reference files and input bam
dir=/data

Time=/usr/bin/time
# Set the reference fasta
#phase2 reference
ref=hs37d5.fa

# Create Variable for input file
INPUT1=kayCombinedNA12878.bam

# sort reads in picard
$Time java $RAM \
    -jar ~/picard.jar \
    SortSam \
    INPUT=${dir}/$INPUT1 \
    OUTPUT=${dir}/$INPUT1.sorted.bam \
    SORT_ORDER=coordinate \
    > sortedKayCombined.report 2>&1
