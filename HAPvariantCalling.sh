# GATK Haplotype VARIANT CALLING

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

# Create Variable for input file
#INPUT1=$INPUT1.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906
INPUT1=NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam

# For running on multipe samples (a cohort) see:
#http://gatkforums.broadinstitute.org/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode

#snps & indels together
time java $RAM \
-jar ~/GenomeAnalysisTK.jar \
-nct $THREADS \
-R ${dir}/${ref} \
-T HaplotypeCaller \
--genotyping_mode DISCOVERY \
--output_mode EMIT_VARIANTS_ONLY \
-I ${dir}/$INPUT1.bqsr.bam \
-o ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
-stand_emit_conf 10.0 \
-stand_call_conf 30.0 \
> HapVars.report 2>&1
	#-I $INPUT2.recal.bam \
	#-I $INPUT3.recal.bam \
	#-I $INPUT4.recal.bam \
	#-metrics snps.metrics \
	#-L /b37/target_intervals.bed \ intervals to operate over. See docstring

## GATK VARIANT QUALITY SCORE RECALIBRATION
#SNP
java $RAM \
  -Djava.io.tmpdir=/tmp \
  -jar ~/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${dir}/${ref} \
  -input ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
  -nt $THREADS \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${dir}/hapmap_3.3.b37.vcf \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 ${dir}/1000G_omni2.5.b37.vcf \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dir}/dbsnp_137.b37.vcf \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${dir}/1000G_phase1.indels.b37.vcf \
  -an QD \
  -an DP \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode SNP \
  -minNumBad 1000 \
  -recalFile ${dir}/$INPUT.SNP.recal \
  -tranchesFile ${dir}/$INPUT.SNP.tranches \
  -rscriptFile ${dir}/$INPUT.SNP.plots.R \
  > VariantRecalibrator_SNP.report 2>&1
  #-percentBad 0.01 \ depreciated

#Apply Snp Recalibration
java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
  -o ${dir}/$INPUT1.SNP.vqsr.SNP.vcf \
  -R ${dir}/${ref} \
  -nt $THREADS \
  -ts_filter_level 99.0 \
  -tranchesFile ${dir}/$INPUT.SNP.tranches \
  -recalFile ${dir}/$INPUT.SNP.recal \
  -mode SNP \
  > ApplyRecalibration_SNP.report 2>&1

#Indel Recalibration
java $RAM -Djava.io.tmpdir=/tmp \
  -jar ~/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${dir}/${ref} \
  -input ${dir}/$INPUT1.unified.raw.SNP.gatk.vcf \
  -nt $THREADS \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 ${dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
  -an DP \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode INDEL \
  -minNumBad 1000 \
  -recalFile ${dir}/$INPUT.INDEL.recal \
  -tranchesFile ${dir}/$INPUT.INDEL.tranches \
  -rscriptFile ${dir}/$INPUT.INDEL.plots.R \
  > VariantRecalibrator_INDEL.report 2>&1
  #-percentBad 0.01 \depreciated

#Apply Indel Recalibration
java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input ${dir}/$INPUT1.unified.raw.INDEL.gatk.vcf \
  -o ${dir}/$INPUT.vqsr_INDEL.vcf \
  -R ${dir}/${ref} \
  -nt $THREADS \
  -ts_filter_level 99.0 \
  -tranchesFile ${dir}/$INPUT.INDEL.tranches \
  -recalFile ${dir}/$INPUT.INDEL.recal \
  -mode INDEL \
  > ApplyRecalibration_INDEL.report 2>&1