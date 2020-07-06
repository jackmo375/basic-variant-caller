#
#	Testing GSNAP aligner
#
#	http://research-pub.gene.com/gmap/
#
#	gsnap is a *genome alignment* program.
#	gmap is the overall package, and gmap is
#	for mRNA alignment, but gsnap uses part of 
#	this library. 
#
#################################################

#!/bin/bash

# include config files
source ../configs/parameters.cfg

main()
{
	species=lambda-virus

	# if the right index does not exist, build it:
	if [[ ! -d "${dat_dir}/indices/${species}.gsnap" ]]; then
		build_index
	fi

	# then chose and run a pipeline:
	#single_read_align
	paired_end_align

}


#
#  functions
#
build_index()
{
	echo 'building gsnap index...'
	${gsnap}/gmap_build \
		-D ${dat_dir}/indices \
		-d ${species}.gsnap \
		${dat_dir}/references/${species}.fa
}

single_read_align()
{
	# Usage: gsnap [OPTIONS...] <FASTA file>
	#
	# option -A specifies the output format

	echo 'aligning single read file...'
	${gsnap}/gsnap \
		-D ${dat_dir}/indices \
		-d ${species}.gsnap \
		-A sam \
		${dat_dir}/reads/${species}_1.fq \
		> ${sam_dir}/${species}.gns.sr.sam
}

paired_end_align()
{
	# So I think, as with most of these tools,
	# if you just specify two input files rather
	# than one, gsnap assumes they are paired
	# end reads...

	echo 'aligning single read file...'
	${gsnap}/gsnap \
		-D ${dat_dir}/indices \
		-d ${species}.gsnap \
		-A sam \
		${dat_dir}/reads/${species}_1.fq \
		${dat_dir}/reads/${species}_2.fq \
		> ${sam_dir}/${species}.gns.pe.sam
}

#
#  run main
#
main
