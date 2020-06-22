#
#	Testing the BWA-MEM aligner
#
#	http://bio-bwa.sourceforge.net/
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
	build_index

	# then run a pipeline
	#simple_align

}


#
#  functions
#
build_index()
{
	echo 'building bwa index...'
	${bwa}/bwa index \
		-p ${dat_dir}/indices/${species}.bwa \
		${dat_dir}/references/${species}.fa
}

simple_align()
{
	echo 'simple align with bwa'

	${bwa}/bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz

}

#
#  run main
#
main
