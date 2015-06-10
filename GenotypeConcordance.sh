#! /bin/bash

#This script uses the gatk genotypeConcordance tool to compare to vcf for concordance

set -e
set -x
set -o pipefail

Time=/usr/bin/time
ref=hs37d5.fa
dir=/data
RAM=-Xmx200g
test_set.vcf=${dir}/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam.unified.raw.BOTH.gatk.vcf
truth_set.vcf=

# get/install samtools
sudo apt-get -qy install zlib1g-dev
sudo apt-get -qy install samtools

# get GATK.jar
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/GenomeAnalysisTK.jar

cd ${dir}
# get the truth_set
wget https://s3-us-west-2.amazonaws.com/bd2k-test-data/platinum-genomes-vcf-NA12878_S1.genome.vcf

#phase2 reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa

cd ~/
#Compare the 2 VCFs
$Time java $RAM \
   -jar ~/GenomeAnalysisTK.jar \
   -T GenotypeConcordance \
   -R ${dir}/${ref} \
   --eval ${test_set.vcf} \
   --comp ${truth_set.vcf} \
   -o ConcordanceOutput.grp \
   > Concordance.report 2>&1