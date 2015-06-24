#! /bin/bash

#This script uses the gatk genotypeConcordance tool to compare to vcf for concordance
#Script assumes you have a reference genome and its index file in $dir

set -e
set -x
set -o pipefail

Time=/usr/bin/time
ref=human_g1k_v37.fasta
dir=/data
RAM=-Xmx200g
test_set_vcf=${dir}/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam.unified.raw.BOTH.gatk.vcf
truth_set_vcf=${dir}/NIST_RTG_PlatGen_merged_highconfidence_v0.2_Allannotate.vcf
#old truth: no worky: platinum-genomes-vcf-NA12878_S1.genome.vcf

# get/install samtools
sudo apt-get -qy install zlib1g-dev
sudo apt-get -qy install samtools

# get GATK.jar
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/GenomeAnalysisTK.jar

#Skipping next 3 lines, have put truth-set on /data on podk-1-1
#cd ${dir}
# get the truth_set
#wget https://s3-us-west-2.amazonaws.com/bd2k-test-data/${truth_set_vcf}

cd ~/
#Compare the 2 VCFs
$Time java $RAM \
   -jar ~/GenomeAnalysisTK.jar \
   -T GenotypeConcordance \
   -R ${dir}/${ref} \
   --eval ${test_set_vcf} \
   --comp ${truth_set_vcf} \
   -o ConcordanceOutput.grp \
   > Concordance.report 2>&1