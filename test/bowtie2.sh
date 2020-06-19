#
# testing the bowtie 2 aligner
#
##################################

#!/bin/bash

# include config files
source ../configs/parameters.cfg

main()
{
	species=lambda-virus

	# if the right index does not exist, build it:
	if [[ ! -f "${dat_dir}/indices/${species}.1.bt2" ]]; then
		echo 'build index...'
		build_index
	fi

	## then chose and run a pipeline:
	simple_align
	#paired_end_align
	#local_align
}


#
#  functions
#
# build reference index
build_index()
{
	${bowtie2}/bowtie2-build \
		${dat_dir}/references/${species}.fa \
		${dat_dir}/indices/${species}
}

## align single reads file
simple_align()
{
	${bowtie2}/bowtie2 \
		-x ${dat_dir}/indices/${species} \
		-U ${dat_dir}/reads/${species}_1.fq \
		-S ${sam_dir}/${species}.sam
}

## align paired-end reads
paired_end_align()
{
	${bowtie2}/bowtie2 \
		-x ${dat_dir}/indices/${species} \
		-1 ${dat_dir}/reads/${species}_1.fq \
		-2 ${dat_dir}/reads/${species}_2.fq \
		-S ${sam_dir}/${species}_paired-end.sam
}

## local alignment
local_align()
{
	${bowtie2}/bowtie2 --local \
		-x ${dat_dir}/indices/${species} \
		-U ${dat_dir}/reads/${species}_long.fq \
		-S ${sam_dir}/${species}_local.sam
}


#
#  run main
#
main

