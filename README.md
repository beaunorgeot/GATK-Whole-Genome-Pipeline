#Purpose
The scripts within this repo were created to run the GATK germline variant calling best pratices pipeline. 
The pipeline takes unprocessed bam file(s) as input and outputs the biological genomic variance of the input to a given reference

##Overview

###Setup & Pipeline
There are currently 4 scripts: 

1. **Setup** installs all of the necessary tools and all of the references.
2.  **processingGATK** runs the steps necessary to transforms reads into being ready for variant calling. Flagstat,Sort, Mark Duplicates, Indel Realignment, BQSR
3.  **UGvariantCalling.sh** Unified Genotyper variant calling plus VQSR
4.  **HAPvariantCalling.sh** Haplotyper variant calling plus VQSR

There is no clear answer about which genotyper to use, you need to pick what you depending on your goals. The HC is more accurate for both SNPs and INDELs but is horribly slow and doesn't scale well to multiple samples.

#
##Reference Files
All reference files have come from the GATK Resource bundle which can be found at [B37 Resource Bundle](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37)

* **Reference Genome (GRCh37):**
* human_g1k_v37_decoy.fasta
* **VCF files for RealignerTargetCreator knowns and dbsnp **for BaseRecalibrator.
* known_indel: /data/GATK_Bundle/*Mills_and_1000G_gold_standard.indels.b37.vcf
*known_indel: /data/GATK_Bundle/*1000G_phase1.indels.b37.vcf
*known_dbsnp: /data/GATK_Bundle/dbsnp_137.b37.vcf
* **1000Genomes Background BAMS ran concurrently with UnifiedGenotyper**
CEU
GBR
* **Resource files for VariantRecalibrator_SNP**
*hapmap_3.3.b37.vcf
*1000G_omni2.5.b37.vcf
*1000G_phase1.snps.high_confidence.b37.vcf
* **Resource files for VariantRecalibrator_INDEL**
*Mills_and_1000G_gold_standard.indels.b37.vcf
*1000G_phase1.indels.b37.vcf

###Workflow
The pre-processing stages of Flagstat,Sort,MarkD,Indel,BQSR are to be run on 1 genome at a time
with the output from one stage being piped in as input for the next (w/the exception of Flagstat)
-Variant calling stages, including VQSR, are to be run (normally) on multiple genomes at once
The input for variant calling stages is the final BQSR output. Calling indel and snp variants are independent
processes that do not depend on eachother in any way. They take the same input. This is also true of VQSR

**Steps**

1. Launch a trusty-ubuntu image
1. Install all scripts and make them executable

		chmod u+x <each Script here>

2. If downloading any documents from s3, Create environment variables for your s3 credentials

        export SECRET_KEY=yourSecretKey
        export ACCESS_KEY=yourAccessKey
3. Run GATKsetup.sh

		./GATKsetup.sh
3
. Use your favorite text editor to assign your input file(s) to the INPUT variable in the GATKpipe.sh script

5. At top of the **processing**, **UGcaller/HCcaller** scripts, set the RAM, THREADS, dir, ref, & INPUT variables to reflect your run. The RAM variable should be set to reflect 90% of your machine's RAM. The THREADS variable should reflect the number of cores
4. Run processingGATK.sh

		./processingGATK.sh
		
5. Run either the HAPvariantCalling.sh OR UGvariantCalling.sh scripts

##Tools
* **Process:**
* Flagstat -get simple metrics on the reads
* Sort -sort reads by reference position. Input is bag of reads that have been mapped
* MarkDuplicates -label/remove reads that map to the same reference position, keeping the read w/the highest raw score
* IndelRealignment -locally realign the reads to minimize the number of mismatches to the reference
* BQSR -rescore the bases based on statistical adherence to known snps and indels
* VariantCalling -Search for (ideally) biological differences differences between input genome and reference genome
* Variant Quality Score Recalibration (VQSR) - The GATK callers are designed to be very lenient in calling variants. VQSR is a two step process that uses machine-learning methods to assign a well calibrated probability to each variant call in a raw call set.

###Flags
* **java processes:** 
* Xmx==max memory to assign for task
* **Indel/BQSR:** 
* -T task/process
* -R Reference file
* -I Input 
* -o Output 
* -nt Number of Threads (IndelRealigner)
* -nct Number of CPU threads to allocate per data thread. (BQSR)
* -known Locations in the genome known to contain indels
* -knownSites Locations in the genome known to contain snps
* **UnifiedGenotyper:**
* -L intervals to operate over. Allows you to call variants on only a subset of the genome. Requires a .bed file to point to the desired regions
* -glm Genotype Liklhoods Model. Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together
* -dcov Target coverage threshold for downsampling to coverage. downsample_to_coverage controls the maximum depth of coverage at each locus.
[50 for 4x, 200 for >30x WGS or Whole exome]
* -stand_emit_conf & stand_call_conf 
stand_emit_conf 10.0 means that it won’t report any potential SNPs with a quality below 10.0; 
but unless they meet the quality threshold set by -stand_call_conf (30.0, in this case), 
they will be listed as failing the quality filter. My choices are 'generic'. 
For variant optimization see: http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk


