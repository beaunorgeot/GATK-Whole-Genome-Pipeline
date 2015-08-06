# ADAM-HapCaller Pipeline

###Purpose
Preprocess bam files using adam, call variants w/GATK's Haplotype caller

### Don't Judge Me
This process is just a hack that I've put together to show the utility and interoperability of adam. I'm expecting parts of this process to change, I'm especially hoping that how adam gets converted into bam changes so we can remove some steps. That's why everything is compartmentalized

### Workflow
1. Use the **SpeedyADAM.sh**  script to preprocess the bam files on an EC2 cluster. The current script assumes a cluster of r3.4xlarge machines built using cgcloud. Takes a bam as input; converts to ada,runs markdups,indel,bqsr, then converts adam to bam.
2. Migrate data. I use distcp to put the bqsr.bam to s3, then s3cmd to pull the bqsr.bam down to a single beefy machine that we own.
3. Use **reOrderBam.sh** to cat the bam and order the chromosomes/reads (takes a long time)
4. Use **adam-hap.sh** to sort the reads using picard
5. Use **HAPvariantCalling.sh** to run the all variant calling steps. 

### What to expect
ThIS pipeline takes about 82 hours for 234GB NA12878. 3-4hrs(16,r3.4xlarge machines) running adam on ec2, 7hrs distcp upload, 5hrs s3cmd download,12hrs reorder, 19hrs sort, 35hrs Hap-caller (32core machine). This is approx 2x as fast as running the GATK best practices pipeline on a 32core machine.

### Improvements
1. (Hopeful)If sortting and ordering can be preserved when converting from ADAM--> BAM, we could shave 30 hours off the pipeline time.
2. (In progress) New spark-based uploader being developed at UCSC should reduce the upload time to s3 down from hours to minutes

* If 1 & 2 happen, this pipeline will go from being 2x faster than GATK best practices, to 3x faster
