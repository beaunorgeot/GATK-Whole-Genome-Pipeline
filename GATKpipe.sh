#! /bin/bash

# This script was created to run the GATK germline variant calling best pratices pipeline. 
# The pipeline takes unprocessed bam file(s) as input and outputs the biological genomic variance 
# of the input to a given reference

# Run the GATKsetup.sh script first. See the associated README

# stderr redirected to descriptively named files.report

cd ~

# Set RAM to use for all java processes
RAM=-Xmx200g
# Set number of threads to the number of cores/machine for speed optimization
THREADS=32

# index bam file
samtools index /ephemeral/mnt/NA12878.bam

#START PIPELINE
# calculate flag stats using samtools
time samtools flagstat \
    /ephemeral/mnt/NA12878.bam

# sort reads in picard
time java $RAM \
    -jar ~/picard/dist/picard.jar \
    SortSam \
    INPUT=/ephemeral/mnt/NA12878.bam \
    OUTPUT=/ephemeral/mnt/NA12878.sorted.bam \
    SORT_ORDER=coordinate
    > sortReads.report 2>&1
rm -r /ephemeral/mnt/NA12878.bam

# mark duplicates reads in picard
time java $RAM \
    -jar ~/picard/dist/picard.jar \
    MarkDuplicates \
    INPUT=/ephemeral/mnt/NA12878.sorted.bam \
    OUTPUT=/ephemeral/mnt/NA12878.mkdup.bam \
    METRICS_FILE=/dev/null #File to write duplication metrics to, Required
    > markDups.report 2>&1
rm -r /ephemeral/mnt/NA12878.sorted.bam

# GATK Indel Realignment
# There are 2 steps to the realignment process
# 1. Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
# 2. Running the realigner over those intervals (IndelRealigner)

# index reference genome, There's no need for this, we have the .fai file
# samtools faidx /ephemeral/mnt/human_g1k_v37.fasta

# create sequence dictionary for reference genome
# We could skip this step. I've download the b37.dict

java $RAM \
    -jar ~/picard/dist/picard.jar \
    CreateSequenceDictionary \
    REFERENCE=/ephemeral/mnt/human_g1k_v37_decoy.fasta \
    OUTPUT=/ephemeral/mnt/human_g1k_v37_decoy.fasta.dict 

# create target indels (find areas likely in need of realignment)
time java $RAM \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -known /ephemeral/mnt/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -known /ephemeral/mnt/1000G_phase1.indels.b37.vcf \
    -R /ephemeral/mnt/human_g1k_v37_decoy.fasta \
    -I /ephemeral/mnt/NA12878.mkdups.bam \
    -o /ephemeral/mnt/target_intervals.list \
    -nt $THREADS 
    > targetIndels.report 2>&1

# realign reads
time java $RAM \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known /ephemeral/mnt/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -known /ephemeral/mnt/1000G_phase1.indels.b37.vcf \
    -R /ephemeral/mnt/human_g1k_v37_decoy.fasta \
    -targetIntervals /ephemeral/mnt/target_intervals.list \
    -I /ephemeral/mnt/NA12878.mkdup.bam \
    -o /ephemeral/mnt/NA12878.realigned.bam \
    -nt $THREADS 
    > realignIndels.report 2>&1
rm -r /ephemeral/mnt/NA12878.mkdup.bam

# GATK BQSR

# create gatk recalibration table
time java $RAM \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R /ephemeral/mnt/human_g1k_v37_decoy.fasta \
    -I /ephemeral/mnt/NA12878.bam \
    -knownSites /ephemeral/mnt/dbsnp_138.vcf \
    -knownSites /ephemeral/mnt/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -knownSites /ephemeral/mnt/1000G_phase1.indels.b37.vcf \
    -o /ephemeral/mnt/recal_data.table \
    -nct 32
    > recalibrationTable.report 2>&1

# recalibrate reads
time java $RAM \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R /ephemeral/mnt/human_g1k_v37_decoy.fasta \
    -I /ephemeral/mnt/NA12878.realigned.bam \
    -BQSR /ephemeral/mnt/recal_data.table \
    -o /ephemeral/mnt/NA12878.bqsr.bam \
    -nct 32
    > Bqsr.report 2>&1

# GATK VARIANT CALLING
# two stages, check for snps, check for indels
# When these stages are run for multiple bams, uncomment the desired number of -I slots
# 1000Genomes Background BAMS can be run concurrently with UnifiedGenotyper to improve accuracy? (CEU,GBR are the recommended ones)

#snps
time java $RAM
	-jar /GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
	-nt $THREADS \
	-R /b37/human_g1k_v37_decoy.fasta \
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
	> SnpVars.report 2>&1

#indels
time java $RAM
	-jar /GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
	-nt $THREADS \
	-R /b37/human_g1k_v37_decoy.fasta \
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
	> indelVars.report 2>&1

## GATK VARIANT QUALITY SCORE RECALIBRATION
# Snp Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T VariantRecalibrator
  -R human_g1k_v37_decoy.fasta
  -input Sample1.raw.vcf
  -nt $THREADS
  -resource: hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf
  -resource: omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.sites.vcf
  -resource: dbsnp,known=false,training=true,truth=false,prior=6.0 dbsnp_135.b37.vcf
  -use_annotation: QD
  -use_annotation: HaplotypeScore
  -use_annotation: MQRankSum
  -use_annotation: ReadPosRankSum
  -use_annotation: FS
  -numBadVariants: 5000
  -mode SNP                           
  -recalFile Sample1_SNP.recal            
  -tranchesFile Sample1_SNP.tranches      
  -rscriptFile Sample1_SNP.plots.R       
  > VariantRecalibrator_SNP.report 2>&1

#Apply Snp Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T ApplyRecalibration
  -input Sample1_raw_SNP.vcf
  -o Sample1_SNP_vqsr_SNP.vcf
  -R human_g1k_v37_decoy.fasta
  -nt $THREADS
  -ts_filter_level: 99.0
  -excludeFiltered : TRUE
  -tranchesFile Sample1.tranches
  -recalFile Sample1.recal
  -mode SNP 
  > ApplyRecalibration_SNP.report 2>&1

#Indel Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T VariantRecalibrator
  -R human_g1k_v37_decoy.fasta
  -input Sample1.raw.vcf
  -nt $THREADS
  -resource: mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf 
  -resource: 1000G,known=false,training=true,truth=true,prior=10.0 1000G_phase1.indels.b37.vcf
  -use_annotation: MQRankSum
  -use_annotation: ReadPosRankSum
  -use_annotation: FS
  -numBadVariants: 5000
  -mode INDEL
  -recalFile Sample1_INDEL.recal            
  -tranchesFile Sample1_INDEL.tranches      
  -rscriptFile Sample1_INDEL.plots.R       
  > VariantRecalibrator_INDEL.report 2>&1

#Apply Indel Recalibration
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T ApplyRecalibration
  -input Sample1_raw_INDEL.vcf
  -o Sample1_vqsr_INDEL.vcf
  -R human_g1k_v37_decoy.fasta
  -nt $THREADS
  -ts_filter_level: 99.0
  -excludeFiltered : TRUE
  -tranchesFile Sample1.tranches
  -recalFile Sample1.recal
  -mode INDEL 
  > ApplyRecalibration_INDEL.report 2>&1

#SELECT VARIANTS
#These steps remove the background files and output SNP and INDEL files, then combine them into a single VCF file

#Select Snp
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T SelectVariants
  -R human_g1k_v37_decoy.fasta
  --variant SNP_variant 
  -o output_selected_SNP.file 
  > SelectVariants_SNP.report 2>&1

#Select Indel
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T SelectVariants
  -R human_g1k_v37_decoy.fasta
  --variant INDEL_variant 
  -o output_Selected_INDEL.file 
  > SelectVariants_INDEL.report 2>&1

#Combine Variants
java $RAM -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
  -T CombineVariants
  -R human_g1k_v37_decoy.fasta
  --variant INDEL_variant 
  --variant  output_Selected_INDEL.file
  --variant output_Selected_SNP.file
  -o Final_combined_variants.vcf 
  > CombineVarients.report 2>&1
