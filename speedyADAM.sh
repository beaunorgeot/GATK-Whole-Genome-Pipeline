#! /bin/bash

#Purpose: Run the preproccessing genomic stages on ADAM as fast as possible. Hand variant-ready reads over to an external Variant Caller
# consider var for spark home

set -e
set -x
set -o pipefail

executor_memory=110G
export ADAM_DRIVER_MEMORY=${executor_memory}
export ADAM_EXECUTOR_MEMORY=${executor_memory}
export SPARK_HOME="/opt/sparkbox/spark"
export ADAM_OPTS="--conf spark.eventLog.enabled=true --conf spark.worker.timeout=500 --conf spark.driver.maxResultSize=10g"
#export ADAM_OPTS="--conf spark.shuffle.service.enable=true"

Time=/usr/bin/time
export ADAM_HOME=~/adam-distribution_2.10-0.17.0
export hdfs_root=hdfs://spark-master:8020
Input1=FEEDME
#NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906

#Get spark -s3 downloader
curl -O https://s3-us-west-2.amazonaws.com/bd2k-artifacts/adam/spark-s3-downloader-0.3-SNAPSHOT.jar
#Get the ADAM jar
#curl -O https://s3-us-west-2.amazonaws.com/bd2k-artifacts/adam-distribution_2.11-0.17.0-bin.tar.gz
curl -LOv http://search.maven.org/remotecontent?filepath=org/bdgenomics/adam/adam-distribution_2.10/0.17.0/adam-distribution_2.10-0.17.0-bin.tar.gz
tar -xvzf adam-distribution_2.10-0.17.0-bin.tar.gz

# pull $Input1 from s3 using s3-downloader
/opt/sparkbox/spark/bin/spark-submit --executor-memory ${executor_memory} \
    ~/spark-s3-downloader-0.3-SNAPSHOT.jar \
    s3://bd2k-test-data/${Input1}.bam \
    ${hdfs_root}/${Input1}.bam \
    --concat #files must be explicitly cat'd

# download dbsnp file for bqsr step
curl -O https://s3-us-west-2.amazonaws.com/bd2k-test-data/dbsnp132_20101103.vcf.gz
gunzip dbsnp132_20101103.vcf.gz
mv dbsnp132_20101103.vcf dbsnp_132.vcf
#put file into hdfs
/opt/sparkbox/hadoop/bin/hadoop fs -put dbsnp_132.vcf ${hdfs_root}/dbsnp_132.vcf

#preemptively remove pre-existing files
hadoop fs -rm -R \
    ${hdfs_root}/dbsnp_132.var.adam || true

# convert known snps file to adam variants file
$Time ${ADAM_HOME}/bin/adam-submit vcf2adam \
    ${hdfs_root}/dbsnp_132.vcf \
    ${hdfs_root}/dbsnp_132.var.adam \
    -onlyvariants \
    > convert_known_snps_2adam.report 2>&1

#remove known snps vcf
/opt/sparkbox/hadoop/bin/hadoop fs -rm -R \
    ${hdfs_root}/dbsnp_132.vcf

# START PIPELINE

#preemptively remove pre-existing files
hadoop fs -rm -R ${hdfs_root}/$Input1.adam || true

# convert to adam, and remove bam
$Time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/$Input1.bam \
    ${hdfs_root}/$Input1.adam \
    > convert_bam2adam.report 2>&1

/opt/sparkbox/hadoop/bin/hadoop fs -rm -R \
    ${hdfs_root}/$Input1.bam

# sort the file
$Time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/user/${USER}/NA12878.adam \
    ${hdfs_root}/user/${USER}/NA12878.sort.adam \
    -sort_reads \
    > sort.report 2>&1

#remove .adam input file
/opt/sparkbox/hadoop/bin/hadoop fs -rm -R \
    ${hdfs_root}/$Input1.adam

# mark duplicate reads
$Time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/$Input1.sort.adam \
    ${hdfs_root}/$Input1.mkdup.adam \
    -mark_duplicate_reads \
    > markDups.report 2>&1

#remove sort.adam input file
/opt/sparkbox/hadoop/bin/hadoop fs -rm -R \
    ${hdfs_root}/$Input1.sort.adam

# realign indels
$Time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/$Input1.mkdup.adam \
    ${hdfs_root}/$Input1.ri.adam \
    -realign_indels \
    > realignIndels.report 2>&1

#remove mkdup.adam input
/opt/sparkbox/hadoop/bin/hadoop fs -rm -R \
    ${hdfs_root}/$Input1.mkdup.adam

# recalibrate quality scores
$Time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/$Input1.ri.adam \
    ${hdfs_root}/$Input1.bqsr.adam \
    -recalibrate_base_qualities \
    -known_snps ${hdfs_root}/dbsnp_132.var.adam \
    > BQSR.report 2>&1

/opt/sparkbox/hadoop/bin/hadoop fs -rm -R \
    ${hdfs_root}/dbsnp_132.var.adam

# Convert .adam file back to bam for external Genotyper
$Time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/$Input1.bqsr.adam \
    ${hdfs_root}/$Input1.bqsr.bam \
    > convert_adam2bam.report 2>&1

#HERE BEGINS VARIANT CALLING STEPS USING .bqsr.bam as the input

# upload to s3 using spark-loader
$Time /opt/sparkbox/spark/bin/spark-submit --executor-memory ${executor_memory} \
    ~/spark-s3-downloader-0.3-SNAPSHOT.jar \
    ${hdfs_root}/${Input1}.bqsr.bam \
    s3://bd2k-test-data/${Input1}.bqsr.attempt2.bam \
    > upload.report 2>&1

