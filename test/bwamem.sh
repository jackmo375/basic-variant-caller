#
#	Testing the BWA-MEM aligner
#
#	http://bio-bwa.sourceforge.net/
#
#	for more help, see:
#	1. http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day1/Sequence%20Alignment_July2015_ShamithSamarajiwa.pdf
#
######################################

#!/bin/bash

# config files
source ../configs/parameters.cfg

main()
{
	species=lambda-virus
	threads=1

	# if the right index does not exist, build it:
	if [[ ! -f "${dat_dir}/indices/${species}.bwa.amb" ]]; then
		build_index
	fi

	## then run a pipeline
	#single_read_align
	paired_end_align

}


#
#  functions
#
build_index()
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

paired_end_align()
{
	echo 'aligning paired-end reads with bwa...'

	${bwa}/bwa mem \
		${dat_dir}/indices/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		${dat_dir}/reads/${species}_2.fq \
		> ${sam_dir}/${species}.bwa.pe.sam

}

#
#  run main
#
main

