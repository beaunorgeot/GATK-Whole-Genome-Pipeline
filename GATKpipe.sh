#! /bin/bash

# This script was created to run the GATK germline variant calling best pratices pipeline. 
# The pipeline takes unprocessed bam file(s) as input and outputs the biological genomic variance 
# of the input to a given reference

# See the associated setup.sh and README

cd /mnt2
# default 
# change to ubuntu linux ami /ephemeral/mnt

# get NA12878 from 1000g from BD2K's s3
~/s3cmd/s3cmd get \
    --secret_key=${SECRET_KEY} \
    --access_key=${ACCESS_KEY} \
    s3://bd2k-test-data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam
mv NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam NA12878.bam

# get reference genome                                                                                                                                                                                        
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/Homo_sapiens_assembly19.fasta.fai
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/Homo_sapiens_assembly19.fasta

# get Known indels
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/1000G_phase1.indels.hg19.sites.fixed.vcf

# get db snp known sites
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/dbsnp_138.hg19.fixed.vcf
mv dbsnp_138.hg19.fixed.vcf dbsnp_138.vcf


cd ~
#TOOLS
# get/install samtools
sudo apt-get install zlib1g-dev
sudo apt-get install samtools

# get/install picard
wget https://github.com/broadinstitute/picard/releases/download/1.130/picard-tools-1.130.zip
unzip picard-tools-1.130.zip

# get GATK.jar
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/GenomeAnalysisTK.jar

# index bam file
samtools index /mnt2/NA12878.bam

#START PIPELINE
# calculate flag stats using samtools
time samtools flagstat \
    /mnt2/NA12878.bam

# sort reads in picard
time java -Xmx200g \
    -jar ~/picard/dist/picard.jar \
    SortSam \
    INPUT=/mnt2/NA12878.bam \
    OUTPUT=/mnt2/NA12878.sorted.bam \
    SORT_ORDER=coordinate
rm -r /mnt2/NA12878.bam

# mark duplicates reads in picard
time java -Xmx200g \
    -jar ~/picard/dist/picard.jar \
    MarkDuplicates \
    INPUT=/mnt2/NA12878.sorted.bam \
    OUTPUT=/mnt2/NA12878.mkdup.bam \
    METRICS_FILE=/dev/null #File to write duplication metrics to, Required
rm -r /mnt2/NA12878.sorted.bam

# GATK Indel Realignment
# There are 2 steps to the realignment process
# 1. Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
# 2. Running the realigner over those intervals (IndelRealigner)

# index reference genome, There's no need for this, we have the .fai file
# samtools faidx /mnt2/human_g1k_v37.fasta

# create sequence dictionary for reference genome
# We could skip this and just download the hg19.dict
# https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/Homo_sapiens_assembly19.dict
java -Xmx200g \
    -jar ~/picard/dist/picard.jar \
    CreateSequenceDictionary \
    REFERENCE=/mnt2/Homo_sapiens_assembly19.fasta \
    OUTPUT=/mnt2/Homo_sapiens_assembly19.fasta.dict 

# create target indels (find areas likely in need of realignment)
time java -Xmx200g \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -known /mnt2/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf \
    -known /mnt2/1000G_phase1.indels.hg19.sites.fixed.vcf \
    -R /mnt2/Homo_sapiens_assembly19.fasta \
    -I /mnt2/NA12878.mkdups.bam \
    -o /mnt2/target_intervals.list \
    -nt 32 #number of threads/cores for task

# realign reads
time java -Xmx200g \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known /mnt2/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf \
    -known /mnt2/1000G_phase1.indels.hg19.sites.fixed.vcf \
    -R /mnt2/Homo_sapiens_assembly19.fasta \
    -targetIntervals /mnt2/target_intervals.list \
    -I /mnt2/NA12878.mkdup.bam \
    -o /mnt2/NA12878.realigned.bam \
    -nt 32 
rm -r /mnt2/NA12878.mkdup.bam

# GATK BQSR

# create gatk recalibration table
time java -Xmx200g \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R /mnt2/Homo_sapiens_assembly19.fasta \
    -I /mnt2/NA12878.bam \
    -knownSites /mnt2/dbsnp_138.vcf \
    -knownSites /mnt2/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf \
    -knownSites /mnt2/1000G_phase1.indels.hg19.sites.fixed.vcf \
    -o /mnt2/recal_data.table \
    -nct 32

# recalibrate reads
time java -Xmx200g \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R /mnt2/Homo_sapiens_assembly19.fasta \
    -I /mnt2/NA12878.realigned.bam \
    -BQSR /mnt2/recal_data.table \
    -o /mnt2/NA12878.bqsr.bam \
    -nct 32

# GATK VARIANT CALLING
# two stages, check for snps, check for indels
# When these stages are run for multiple bams, uncomment the desired number of -I slots

#snps
time java -Xmx200g
	-jar /GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
	-nt 32 \
	-R /b37/Homo_sapiens_assembly19.fasta \
	-T UnifiedGenotyper \
	-I file1.recal.bam \
	#-I file2.recal.bam \
	#-I file3.recal.bam \
	#-I file4.recal.bam \
	-o unified.SNP.gatk.vcf \
	#-metrics snps.metrics \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-dcov 200 \
	#-L /b37/target_intervals.bed \ intervals to operate over. See docstring
	-glm SNP \

#indels
time java -Xmx200g
	-jar /GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
	-nt 32 \
	-R /b37/Homo_sapiens_assembly19.fasta \
	-T UnifiedGenotyper \
	-I file1.recal.bam \
	#-I file2.recal.bam \
	#-I file3.recal.bam \
	#-I file4.recal.bam \
	-o unified.INDEL.gatk.vcf \
	#-metrics snps.metrics \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-dcov 200 \
	#-L /b37/target_intervals.bed \ #see docstring
	-glm INDEL \

## NEXT APPROPRIATE STEP WOULD BE VARIANT QUALITY SCORE RECALIBRATION







