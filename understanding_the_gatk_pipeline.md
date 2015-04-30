#Understanding the GATK Pipeline

##Sorting & Mark Duplicates
During the sequencing process, the same DNA molecules can be sequenced several times.The resulting duplicate reads are not informative and should not be counted as additional evidence for or against a putative variant. The duplicate-marking process (sometimes called “dedupping” in bioinformatics slang) identiﬁes these reads as such so that the GATK tools know to ignore them

###If you have a sam input, you can convert to bam, sort, and mark in 1 step

		java -jar MarkDuplicates.jar
		INPUT=aligned reads.sam
		OUTPUT=dedup reads.bam
		SO=coordinate

##INDEL Realignment
The mapping algorithms that are used in the initial step of aligning the data to the reference are prone to various types of artifacts. For example, reads that align on the edges of indels often get mapped with mismatching bases that might look like evidence for SNPs, but are actually mapping artifacts. The realignment process identiﬁes the most consistent placement of the reads with respect to the indel in order to clean up these artifacts

##BQSR
The recalibration process applies an empirically accurate error model to the bases, producing a BAM ﬁlethat is suitable for analysis.
This creates a GATKReport ﬁle called recal data.grp containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data. This also generates a plot called before recal.pdf showing how the original base qualities in your data differ from the base qualities empirically calculated by the Base Recalibrator.

##Haplotype Genotyper (Variant Caller)
The end product of this protocol will be a VCF ﬁle containing raw calls that should be ﬁltered before they can be used in downstream analyses.

Many variant callers specialize in either SNPs or indels, or (like the GATK’s ownUniﬁedGenotyper) have to call them using separate models of variation. 

The Haplotype Caller is capable of calling SNPs and indels simultaneously via local de novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation,it discards the existing mapping information and completely reassembles the reads in that region. This allows the Haplotype Caller to be more accurate when calling regions that are traditionally difﬁcult to call, for example,when they contain different types of variants close to each other. It also makes the HC much better at calling indels.

You can specify how you want the program to determine the alternate alleles to use for genotyping. In the default DISCOVERY mode, the program will choose the mostlikely alleles out of those it sees in the data. If you just want to determine if a sample has a speciﬁc genotype of interest and you are not interested in other alleles, you could change the mode to GENOTYPE GIVEN ALLELES and pass a VCF (w/the -alleles option) to check for the presence of specific alleles.

Emission conﬁdence threshold (--stand emit conf) and Call Confidence threshold (-stand call conf)
Emit threshold sets the lower limit for the score of any base that should even be considered for the possibility of being a variant. Call threshold sets the lower limit for any base that is actually determined to be a variant. 
I THINK these scores actually relate to the HC’s confidence in whether a base is actually a variant, which is not synonymous with the sequencing confidence for a base. Consider the case where there are many reads covering a position, all of them are G’s, but all of them have low sequencing confidence

Although, at the end of this step, you now have a nice fresh set of variant calls, the variant-discovery stage is not over. The distinctions made by the caller itself between low-conﬁdence calls and the rest are very primitive, and should not be taken as a deﬁnitive guide for ﬁltering. The GATK callers are designed to be very lenient in calling variants, so it is extremely important to apply one of the ﬁltering methods

##VQSR
In the ﬁrst step of this two-step process, the program uses machine-learning methods to assign a well calibrated probability to each variant call in a raw call set. We can then use this variant quality score in the second step to ﬁlter the raw call set, thus producing a  subset of calls with our desired level of quality, ﬁne-tuned to balance speciﬁcity and sensitivity.

To calculate the variant quality score, the program builds an adaptive error model using training sets (explained further below). The model provides a continuous, covaryingestimate of the relationship between variant call annotations and the probability that avariant call is a true genetic variant, rather than a sequencing or data-processing artifact.The program then applies this adaptive error model to both known and novel variationdiscoveredinthecallsetofinterest,andannotateseachvariantwithaqualityscorecalledVQSLOD. This is the log odds ratio of the variant call being a true positive versus being a false positive according to the training model.

For each training set,we use key-value tags to qualify whether the set contains known sites, training sites, and/or truth sites. We also use a tag to specify the prior likelihood that those sites are true (using the Phred scale).

Specify which annotations the program should use to evaluate the likelihood of SNPs being real. These annotations are included in the information generated foreach variant call by the caller.

##APPLY RECALIBRATION
This creates a new VCF ﬁle,which contains all the original variants from the original raw variants.vcf ﬁle, but now the SNPs are annotated with their recalibrated quality scores(VQSLOD) and either PASS or  FILTER depending on whether or not they are included in the selected tranche.

It's recommended that you take the second lowest of the tranches speciﬁed in the original recalibration command. This means that we are applying to our data set the level of sensitivity that would allow us to retrieve 99% of true variants from the truth-training sets of HapMap and Omni SNPs.If we wanted to be more speciﬁc(and therefore have less risk of including false positives, at the risk of missing real sites), we could take the very lowest tranche,which would only retrieve 90% of the truth-training sites.If we wanted to be more sensitive(and therefore less speciﬁc, at the risk of including more false positives), we could take the higher tranches. In our Best Practices documentation, we recommend taking the second-highest tranche (99.9%), which provides the highest sensitivity you can get while still being acceptably speciﬁc.


