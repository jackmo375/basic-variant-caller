#
#  testing the full BWA MEM + GATK pipeline
#
#	https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4
#	http://bio-bwa.sourceforge.net/
#
#  formerly, this pipeline is for:
#	*germline short variant discovery (SNPs and Indels)*
#
#	see broad institute for best practices:
#	https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
#
#######################################################

#!/bin/bash

source ../configs/parameters.cfg

main()
{
	species=lambda-virus
	threads=1

	# if the right index does not exist, build it:
	if [[ ! -f "${dat_dir}/indices/${species}.bwa.amb" ]]; then
		build_align_index
	fi

	build_varcall_dict
	#build_varcall_index

	## then run the pipeline steps
	single_read_align
	#call_variants
}

#
#  functions
#
build_align_index()
{
	# option -a specifies the algorithm, with arguments:
	# 
	# 1. 'is' is for short sequences (like virus genomes)
	# 2. 'bwtsw' is for long sequences (like human genomes)
	
	echo 'building bwa index...'
	${bwa}/bwa index \
		-p ${dat_dir}/indices/${species}.bwa \
		-a is \
		${dat_dir}/references/${species}.fa
}

build_varcall_dict()
{

}

single_read_align()
{
	echo 'aligning single reads with bwa...'

	# basic usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	# <idxbase> :: the base of the index file names
	#				i.e. the names *with* full path, but no extension
	# <in1.fq>  :: input file of sequence reads in fastq/fasta format
	# [in2.fq]  :: (optional) input file of mates for the first read file, 
	#				if the input data is paired-end 

	${bwa}/bwa mem \
		${dat_dir}/indices/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		> ${sam_dir}/${species}.bwa.sr.sam
}

call_variants()
{
	echo 'calling variants with gatk...'

	${gatk}/gatk HaplotypeCaller \
		--input ${sam_dir}/${species}.bwa.sr.sam \
		--output ${vcf_dir}/${species}.bwa-gatk.sr.vcf \
		--reference ${dat_dir}/references/${species}.fa
}

#
#  run main
#
main
