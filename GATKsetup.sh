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