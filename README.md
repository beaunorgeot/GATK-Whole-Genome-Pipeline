#Purpose
This script was created to run the GATK germline variant calling best pratices pipeline. 
The pipeline takes unprocessed bam file(s) as input and outputs the biological genomic variance of the input to a given reference

##Overview
###Setup & Pipeline
There are 2 scripts, Setup and Pipeline. The setup script installs all of the necessary tools and all of the references. The pipeline script implements the GATK pipeline.
###Workflow
The pre-processing stages of Flagstat,Sort,MarkD,Indel,BQSR are to be run on 1 genome at a time
with the output from one stage being piped in as input for the next (w/the exception of Flagstat)
-Variant calling stages, including VQSR, are to be run (normally) on multiple genomes at once
The input for variant calling stages is the final BQSR output. Calling indel and snp variants are independent
processes that do not depend on eachother in any way. They take the same input. This is also true of VQSR

##Tools
* **Process:**
* Flagstat -get simple metrics on the reads
* Sort -sort reads by reference position. Input is bag of reads that have been mapped
* MarkDuplicates -label/remove reads that map to the same reference position, keeping the read w/the highest raw score
* IndelRealignment -locally realign the reads to minimize the number of mismatches to the reference
* BQSR -rescore the bases based on statistical adherence to known snps and indels
* VariantCalling -Search for (ideally) biological differences differences between input genome and reference genome

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

##Warnings:
-Make sure to place your actual bam files into -I flag Locations
-Update -L flags with the target_intervals for each bam
-secret_keys
-Change nt/nct appropriately for each instance type. 
