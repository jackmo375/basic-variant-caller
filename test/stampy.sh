#
#	testing stampy aligner
#
#	https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy
#
############################

#!/bin/bash

# include parameter files

source ../configs/parameters.cfg

main()
{
	species=lambda-virus

	if [[ ! -f "${dat_dir}/indices/${species}.sta.stidx" ]]; then
		build_index
	fi

	# then chose a pipeline:
	# (note: paired_end_align is not yet working)
	single_read_align
	#paired_end_align

}


#
#  functions
#
build_index()
{
	#  for Stampy you need to build both an index file (stidx)
	#  and a hash table (.sthash) from a reference genome before
	#  performing any alignment. this function does both together. 
	#
	# the options are a little hard to understand, but something like this
	# uppercase options (-G, -H) specify the file base of the outputs
	# lowercase options (-g -h) specify the base of the inputs.
	#
	# the file base is the name of the file, *without* the extension. 
	# stampy accepts both a relative- or absolute-path file base. 

	# g/G specifies an index file, and h/H is a hash table. So, if we specify
	# that an input is an index file (-g) and an output is a hash table (-H) 
	# for stampy.py then it is clear what stampy.py is meant to do here! Clever. 

	echo 'building stampy index...'

	# first we build an index (.stidx)
	${stampy}/stampy.py \
		--overwrite \
		-G ${dat_dir}/indices/${species}.sta \
		${dat_dir}/references/${species}.fa

	# build hash table
	${stampy}/stampy.py \
		--overwrite \
		-g ${dat_dir}/indices/${species}.sta \
		-H ${dat_dir}/indices/${species}.sta

		
}

single_read_align()
{
	echo 'performing single read alignment...'

	${stampy}/stampy.py \
		-g ${dat_dir}/indices/${species}.sta \
		-h ${dat_dir}/indices/${species}.sta \
		-M ${dat_dir}/reads/${species}_1.fq \
		-f sam \
		> ${sam_dir}/${species}.sta.sr.sam
}

paired_end_align()
{
	echo 'performing paired end alignment...'

	${stampy}/stampy.py \
		-g ${dat_dir}/indices/${species}.sta \
		-h ${dat_dir}/indices/${species}.sta \
		-M ${dat_dir}/reads/${species}_1.fq,${dat_dir}/reads/${species}_2.fq \
		-f sam \
		> ${sam_dir}/${species}.sta.pe.sam
}

#
#  run main
#
main
