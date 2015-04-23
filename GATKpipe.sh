#! /bin/bash

# This script was created to run the GATK germline variant calling best pratices pipeline. 
# The pipeline takes unprocessed bam file(s) as input and outputs the biological genomic variance 
# of the input to a given reference

# Run the GATKsetup.sh script first. See the associated README

# stderr redirected to descriptively named files.report

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
# Set the reference fasta
# CHECK LINE 83 OR THERE-ABOUTS, IT POINTS TO A REFERENCE NOT INCLUDED IN THIS VAR
ref=hs37d5.fa

cd ${dir}
# get input bam file
#~/s3cmd/s3cmd get s3://bd2k-test-data/$INPUT1.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam

#NA12878 chr20 bam
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam 

# Create Variable for input file
#INPUT1=$INPUT1.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906 
INPUT1=NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam

# index bam file
samtools index ${dir}/$INPUT1

cd ~

#START PIPELINE
# calculate flag stats using samtools
time samtools flagstat \
    ${dir}/$INPUT1

# sort reads in picard
time java $RAM \
    -jar ~/picard/picard.jar \
    SortSam \
    INPUT=${dir}/$INPUT1 \
    OUTPUT=${dir}/$INPUT1.sorted.bam \
    SORT_ORDER=coordinate 
    > sortReads.report 2>&1
rm -r ${dir}/$INPUT1

# mark duplicates reads in picard
time java $RAM \
    -jar ~/picard/picard.jar \
    MarkDuplicates \
    INPUT=${dir}/$INPUT1.sorted.bam \
    OUTPUT=${dir}/$INPUT1.mkdups.bam \
    METRICS_FILE=metrics.txt \
    ASSUME_SORTED=true  
    > markDups.report 2>&1
rm -r ${dir}/$INPUT1.sorted.bam

# GATK Indel Realignment
# There are 2 steps to the realignment process
# 1. Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
# 2. Running the realigner over those intervals (IndelRealigner)

# index reference genome, There's no need for this, we have the .fai file
# samtools faidx ${dir}/${ref}

# create sequence dictionary for reference genome
# We could skip this step. I've download the b37.dict

java $RAM \
    -jar ~/picard/picard.jar \
    CreateSequenceDictionary \
    REFERENCE=${dir}/${ref} \
    OUTPUT=${dir}/hs37d5.dict
    #OUTPUT=${dir}/human_g1k_v37.dict 
    #OUTPUT=${dir}/${ref}.dict

# Index the markdups.bam for realignment
samtools index ${dir}/$INPUT1.mkdups.bam

# create target indels (find areas likely in need of realignment)
time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -known ${dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -known ${dir}/1000G_phase1.indels.b37.vcf \
    -R ${dir}/${ref} \
    -I ${dir}/$INPUT1.mkdups.bam \
    -o ${dir}/target_intervals.list \
    -nt $THREADS 
    > targetIndels.report 2>&1

# realign reads
time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known ${dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -known ${dir}/1000G_phase1.indels.b37.vcf \
    -R ${dir}/${ref} \
    -targetIntervals ${dir}/target_intervals.list \
    -I ${dir}/$INPUT1.mkdups.bam \
    -o ${dir}/$INPUT1.realigned.bam \
    #-nt $THREADS \ nt is not supported by IndelRealigner
    > realignIndels.report 2>&1
rm -r ${dir}/$INPUT1.mkdups.bam

# GATK BQSR

# create gatk recalibration table
time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${dir}/${ref} \
    -I ${dir}/$INPUT1 \
    -knownSites ${dir}/dbsnp_137.vcf \
    -knownSites ${dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -knownSites ${dir}/1000G_phase1.indels.b37.vcf \
    -o ${dir}/recal_data.table \
    -nct $THREADS
    > recalibrationTable.report 2>&1

# recalibrate reads
time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${dir}/${ref} \
    -I ${dir}/$INPUT1.realigned.bam \
    -BQSR ${dir}/recal_data.table \
    -o ${dir}/$INPUT1.bqsr.bam \
    -nct $THREADS \
    > Bqsr.report 2>&1

# GATK VARIANT CALLING
# two stages, check for snps, check for indels
# When these stages are run for multiple bams, uncomment the desired number of -I slots
# 1000Genomes Background BAMS can be run concurrently with UnifiedGenotyper to improve accuracy? (CEU,GBR are the recommended ones)

#snps
time java $RAM
	-jar ~/GenomeAnalysisTK.jar \
	-nt $THREADS \
	-R /b37/${ref} \
	-T UnifiedGenotyper \
	-I $INPUT1.bqsr.bam \
	#-I $INPUT2.recal.bam \
	#-I $INPUT3.recal.bam \
	#-I $INPUT4.recal.bam \
	-o $INPUT1.unified.raw.SNP.gatk.vcf \
	#-metrics snps.metrics \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-dcov 200 \
	#-L /b37/target_intervals.bed \ intervals to operate over. See docstring
	-glm SNP \
	> SnpVars.report 2>&1

#indels
time java $RAM
	-jar ~/GenomeAnalysisTK.jar \
	-nt $THREADS \
	-R /b37/${ref} \
	-T UnifiedGenotyper \
	-I $INPUT1.bqsr.bam \
	#-I $INPUT2.recal.bam \
	#-I $INPUT3.recal.bam \
	#-I $INPUT4.recal.bam \
	-o $INPUT1.unified.raw.INDEL.gatk.vcf \
	#-metrics snps.metrics \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-dcov 200 \
	#-L /b37/target_intervals.bed \ #see docstring
	-glm INDEL \
	> indelVars.report 2>&1

## GATK VARIANT QUALITY SCORE RECALIBRATION
# Snp Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${ref} \
  -input $INPUT1.unified.raw.SNP.gatk.vcf \
  -nt $THREADS \
  -resource: hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
  -resource: omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.sites.vcf \
  -resource: dbsnp,known=false,training=true,truth=false,prior=6.0 dbsnp_135.b37.vcf \
  -use_annotation: QD \
  -use_annotation: HaplotypeScore \
  -use_annotation: MQRankSum \
  -use_annotation: ReadPosRankSum \
  -use_annotation: FS \
  -numBadVariants: 5000 \
  -mode SNP \                           
  -recalFile $INPUT1_SNP.recal \       
  -tranchesFile $INPUT1_SNP.tranches \     
  -rscriptFile $INPUT1_SNP.plots.R       
  > VariantRecalibrator_SNP.report 2>&1

#Apply Snp Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input $INPUT1.unified.raw.SNP.gatk.vcf \
  -o $INPUT1.SNP.vqsr.SNP.vcf \
  -R ${ref} \
  -nt $THREADS \
  -ts_filter_level: 99.0 \
  -excludeFiltered : TRUE \
  -tranchesFile $INPUT1.tranches \
  -recalFile $INPUT1.recal \
  -mode SNP \
  > ApplyRecalibration_SNP.report 2>&1

#Indel Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${ref} \
  -input $INPUT1.unified.raw.INDEL.gatk.vcf \
  -nt $THREADS \
  -resource: mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf \
  -resource: 1000G,known=false,training=true,truth=true,prior=10.0 1000G_phase1.indels.b37.vcf \
  -use_annotation: MQRankSum \
  -use_annotation: ReadPosRankSum \
  -use_annotation: FS \
  -numBadVariants: 5000 \
  -mode INDEL \
  -recalFile $INPUT1_INDEL.recal \            
  -tranchesFile $INPUT1_INDEL.tranches \      
  -rscriptFile $INPUT1_INDEL.plots.R \       
  > VariantRecalibrator_INDEL.report 2>&1

#Apply Indel Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input $INPUT1.unified.raw.INDEL.gatk.vcf \
  -o $INPUT1_vqsr_INDEL.vcf \
  -R ${ref} \
  -nt $THREADS \
  -ts_filter_level: 99.0 \
  -excludeFiltered : TRUE \
  -tranchesFile $INPUT1.tranches \
  -recalFile $INPUT1.recal \
  -mode INDEL \
  > ApplyRecalibration_INDEL.report 2>&1

#SELECT VARIANTS
#These steps remove the background files and output SNP and INDEL files, then combine them into a single VCF file

#Select Snp
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${ref} \
  --variant SNP_variant \
  -o output_selected_SNP.file \
  > SelectVariants_SNP.report 2>&1

#Select Indel
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${ref} \
  --variant INDEL_variant \
  -o output_Selected_INDEL.file \
  > SelectVariants_INDEL.report 2>&1

#Combine Variants
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar \
  -T CombineVariants \
  -R ${ref} \
  --variant INDEL_variant \
  --variant  output_Selected_INDEL.file \
  --variant output_Selected_SNP.file \
  -o Final_combined_variants.vcf \
  > CombineVarients.report 2>&1
