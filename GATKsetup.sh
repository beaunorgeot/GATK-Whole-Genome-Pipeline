# -e for exit, -x for echo, -o pipefail for report failure if ANY process fails (default is to report success if first)
set -e
set -x 
set -o pipefail

cd /mnt/ephemeral

#get reference genome
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/human_g1k_v37.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/human_g1k_v37.fasta.fai.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/human_g1k_v37.dict.gz
gunzip human_g1k_v37.fasta
gunzip human_g1k_v37.fasta.fai 
gunzip human_g1k_v37.dict.gz

# VCF files for VariantRecalibrator_INDEL.RealignerTargetCreator knowns, and dbsnp for BaseRecalibraton
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/1000G_phase1.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/dbsnp_137.b37.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.b37.vcf.gz
gunzip 1000G_phase1.indels.b37.vcf.gz
gunzip dbsnp_137.b37.vcf.gz

#Resource files for VariantRecalibrator_SNP
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/hapmap_3.3.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/1000G_omni2.5.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
gunzip hapmap_3.3.b37.vcf.gz
gunzip 1000G_omni2.5.b37.vcf.gz
gunzip 1000G_phase1.snps.high_confidence.b37.vcf.gz


cd ~
#TOOLS

# install s3cmd
sudo apt-get -qy install s3cmd

# get/install samtools
sudo apt-get -qy install zlib1g-dev
sudo apt-get -qy install samtools

# get/install picard, requires unzip program
wget https://github.com/broadinstitute/picard/releases/download/1.130/picard-tools-1.130.zip
sudo apt-get install unzip
unzip picard-tools-1.130.zip

# get GATK.jar
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/GenomeAnalysisTK.jar

cat > ~/.s3cfg <<-END
	[default]
	access_key = $ACCESS_KEY
	bucket_location = US
	cloudfront_host = cloudfront.amazonaws.com
	default_mime_type = binary/octet-stream
	delete_removed = False
	dry_run = False
	enable_multipart = True
	encoding = UTF-8
	encrypt = False
	follow_symlinks = False
	force = False
	get_continue = False
	gpg_command = /usr/bin/gpg
	gpg_decrypt = %(gpg_command)s -d --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
	gpg_encrypt = %(gpg_command)s -c --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
	gpg_passphrase =
	guess_mime_type = True
	host_base = s3.amazonaws.com
	host_bucket = %(bucket)s.s3.amazonaws.com
	human_readable_sizes = False
	invalidate_on_cf = False
	list_md5 = False
	log_target_prefix =
	mime_type =
	multipart_chunk_size_mb = 15
	preserve_attrs = True
	progress_meter = True
	proxy_host =
	proxy_port = 0
	recursive = False
	recv_chunk = 4096
	reduced_redundancy = False
	secret_key = $SECRET_KEY
	send_chunk = 4096
	simpledb_host = sdb.amazonaws.com
	skip_existing = False
	socket_timeout = 300
	urlencoding_mode = normal
	use_https = False
	verbosity = WARNING
	website_endpoint = http://%(bucket)s.s3-website-%(location)s.amazonaws.com/
	website_error =
	website_index = index.html
END