#
#  gatk
#
#	simple pipeline for
#		+ simulated SNP & indel data
#		+ short virus genomes
#
#	Jack Morrice
#
#	we are using this workflow for inspiration:
#	https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.wdl
#
###############################################

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { argv=("$@")
	
	pipeline=gatk
	input_json=${argv[0]}
	output_label=$(basename $input_json .json)

	# get pipeline arguments from input pipeline json

	species=$(value_from_json $input_json '.species')
	threads=$(value_from_json $input_json '.threads')
	read_type=$(value_from_json $input_json '.read_type')

	# execute workflow steps

	if [[ ${read_type} == 'pe' ]]; then 
		custom_call paired_end_align "aligning paired-end reads with bwa..."
	else
		custom_call single_read_align "aligning reads to the reference..."
	fi

	custom_call sam_to_bam "converting sam files to bam format..."

	custom_call sort_sam "sorting bam files..."

	custom_call mark_duplicates "marking duplicates..."

	custom_call validate_bam "validating bam files..."

	## for now we will skip the recalibration of base quality scores. 

	custom_call index_bam "indexing all input bam files..."

	custom_call call_variants "calling variants..."

}	


#
#  tasks
#
single_read_align()
{
	# basic usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	# <idxbase> :: the base of the index file names
	#				i.e. the names *with* full path, but no extension
	# <in1.fq>  :: input file of sequence reads in fastq/fasta format
	# [in2.fq]  :: (optional) input file of mates for the first read file, 
	#				if the input data is paired-end 
	#
	# gatk requires read group information. I am not sure exactly what
	# this means, but there is information here:
	# https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups

	$bwa mem \
		${ref_dir}/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		-R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:LIB-DAD-1\tSM:DAD\tPI:200" \
		-t $threads \
		> ${sam_dir}/${output_label}.sam
}

paired_end_align()
{
	${bwa} mem \
		${ref_dir}/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		${dat_dir}/reads/${species}_2.fq \
		-R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:LIB-DAD-1\tSM:DAD\tPI:200" \
		-t $threads \
		> ${sam_dir}/${output_label}.sam

}

sam_to_bam()
{
	${samtools} view \
		-S \
		-b ${sam_dir}/${output_label}.sam \
		> ${bam_dir}/${output_label}.bam
}

sort_sam()
{
	# the sam files need to be sorted before we can mark any duplicates. 
	# this can be done with either picard or samtools.
	#
	# even though the picard function is called 'SortSam' it seems to 
	# actually sort bam files. That's fine. 

	java -jar ${picard} SortSam \
		I=${bam_dir}/${output_label}.bam \
		O=${bam_dir}/${output_label}.sorted.bam \
		SORT_ORDER=coordinate
}

mark_duplicates()
{
	# this can be done using either samtools or picard.
	# http://www.htslib.org/algorithms/duplicate.html
	# https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

	java -jar ${picard} MarkDuplicates \
		I=${bam_dir}/${output_label}.sorted.bam \
		O=${bam_dir}/${output_label}.marked.bam \
		M=${bam_dir}/${output_label}.marked_dup_metrics.txt

}

validate_bam()
{
	java -jar $picard ValidateSamFile \
		I=${bam_dir}/${output_label}.marked.bam \
		MODE=SUMMARY
}

index_bam()
{
	$samtools index \
		${bam_dir}/${output_label}.marked.bam

}

call_variants()
{
	# I think we need to perform variant calling for a sample ensemble?
	# doesn't seem to work for a single sample? Not sure...

	$gatk --java-options "-Xmx4g" HaplotypeCaller  \
		--native-pair-hmm-threads $threads \
		-R ${ref_dir}/${species}.fa \
		-I ${bam_dir}/${output_label}.marked.bam \
		-O ${vcf_dir}/${output_label}.raw.g.vcf
}



#
#  run workflow
#
workflow "$@"

