#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { argv=("$@")

	species=${argv[0]}	

	custom_call build_index "building alignment index..."

	custom_call build_varcall_index "building variant caller index..."

	custom_call build_varcall_dict "building variant caller dictionary"
}

#
#  tasks
#
build_index()
{
	# option -a specifies the algorithm, with arguments:
	# 
	# 1. 'is' is for short sequences (like virus genomes)
	# 2. 'bwtsw' is for long sequences (like human genomes)


	if [[ ! -f "${ref_dir}/${species}.bwa.amb" ]]; then
		$bwa index \
			-p ${ref_dir}/${species}.bwa \
			-a is \
			${ref_dir}/${species}.fa
	else
		printf "skipping index build as index already exists"; echo
	fi
	
}

build_varcall_index() 
{
	# the gatk haplotype variant caller requires a dictionary and an index 
	# are both built from the reference. this function builds the index 
	# which I think is different from the alignment index?
	#
	# http://www.htslib.org/doc/samtools-faidx.html
	# https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

	if [[ ! -f "${ref_dir}/${species}.fai" ]]; then
		$samtools faidx \
			${ref_dir}/${species}.fa \
			> ${ref_dir}/${species}.fai
	else
		printf "skipping varcall index build as this index already exists"; echo
	fi

}

build_varcall_dict()
{
	# the gatk haplotype variant caller requires a dictionary and an index 
	# are both built from the reference. this function builds the dictionary 
	# which I think is different from the alignment index?
	#
	# https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
	# https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-

	if [[ ! -f "${ref_dir}/${species}.dict" ]]; then
		java -jar $picard CreateSequenceDictionary \
			R=${ref_dir}/${species}.fa \
			O=${ref_dir}/${species}.dict
	else
		printf "skipping varcall dictionary build as this dictionary already exists"; echo
	fi

}

#
#  run workflow
#
workflow "$@"

