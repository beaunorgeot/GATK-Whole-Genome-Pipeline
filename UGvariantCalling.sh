# GATK Unified Genotyper VARIANT CALLING

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
#phase2 reference
ref=hs37d5.fa
#choose how to log time
Time=/usr/bin/time

# Create Variable for input file
#INPUT1=$INPUT1.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906
INPUT1=NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam

# two stages, check for snps, check for indels
# When these stages are run for multiple bams, uncomment the desired number of -I slots
# 1000Genomes Background BAMS can be run concurrently with UnifiedGenotyper to improve accuracy? (CEU,GBR are the recommended ones)

#snps
$Time java $RAM \
-jar ~/GenomeAnalysisTK.jar \
-nt $THREADS \
-R ${dir}/${ref} \
-T UnifiedGenotyper \
-I ${dir}/$INPUT1.bqsr.bam \
-o ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
-stand_call_conf 30.0 \
-stand_emit_conf 10.0 \
-dcov 200 \
-glm SNP \
> SnpVars.report 2>&1
	#-I $INPUT2.recal.bam \
	#-I $INPUT3.recal.bam \
	#-I $INPUT4.recal.bam \
	#-metrics snps.metrics \
	#-L /b37/target_intervals.bed \ intervals to operate over. See docstring

#indels
$Time java $RAM \
-jar ~/GenomeAnalysisTK.jar \
-nt $THREADS \
-R ${dir}/${ref} \
-T UnifiedGenotyper \
-I ${dir}/$INPUT1.bqsr.bam \
-o ${dir}/$INPUT1.unified.raw.INDEL.gatk.vcf \
-stand_call_conf 30.0 \
-stand_emit_conf 10.0 \
-dcov 200 \
-glm INDEL \
> indelVars.report 2>&1
	#-I $INPUT2.recal.bam \
	#-I $INPUT3.recal.bam \
	#-I $INPUT4.recal.bam \
	#-metrics snps.metrics \
	#-L /b37/target_intervals.bed \ #see docstring

## GATK VARIANT QUALITY SCORE RECALIBRATION
# Snp Recalibration
$Time java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${dir}/${ref} \
  -input ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
  -nt $THREADS \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${dir}/hapmap_3.3.b37.vcf \
  -resource:omni,known=false,training=true,truth=false,prior=12.0 ${dir}/1000G_omni2.5.b37.vcf \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dir}/dbsnp_137.b37.vcf \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${dir}/1000G_phase1.indels.b37.vcf \
  -an QD \
  -an DP \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode SNP \
  -recalFile ${dir}/$INPUT1.SNP.recal \
  -tranchesFile ${dir}/$INPUT1.SNP.tranches \
  -rscriptFile ${dir}/$INPUT1.SNP.plots.R \
  > VariantRecalibrator_SNP.report 2>&1
#try simply removing the additional args at the end of 195:197 completely?
# But my current layout is exactly like the docs:https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php


#Apply Snp Recalibration
$Time java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
  -o ${dir}/$INPUT1.SNP.vqsr.SNP.vcf \
  -R ${dir}/${ref} \
  -nt $THREADS \
  -ts_filter_level 99.0 \
  -tranchesFile ${dir}/$INPUT1.SNP.tranches \
  -recalFile ${dir}/$INPUT1.SNP.recal \
  -mode SNP \
  > ApplyRecalibration_SNP.report 2>&1
  # removed:  -excludeFiltered : TRUE \


#Indel Recalibration
$Time java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${dir}/${ref} \
  -input ${dir}/$INPUT1.unified.raw.INDEL.gatk.vcf \
  -nt $THREADS \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 ${dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
  -an DP \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode INDEL \
  -recalFile ${dir}/$INPUT1.INDEL.recal \
  -tranchesFile ${dir}/$INPUT1.INDEL.tranches \
  -rscriptFile ${dir}/$INPUT1.INDEL.plots.R \
  > VariantRecalibrator_INDEL.report 2>&1

#Apply Indel Recalibration
$Time java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input ${dir}/$INPUT1.unified.raw.INDEL.gatk.vcf \
  -o ${dir}/$INPUT1.vqsr_INDEL.vcf \
  -R ${dir}/${ref} \
  -nt $THREADS \
  -ts_filter_level 99.0 \
  -tranchesFile ${dir}/$INPUT1.INDEL.tranches \
  -recalFile ${dir}/$INPUT1.INDEL.recal \
  -mode INDEL \
  > ApplyRecalibration_INDEL.report 2>&1

