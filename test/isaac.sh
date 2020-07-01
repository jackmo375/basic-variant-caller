#
#	testing isaac aligner
#
#	https://github.com/Illumina/Isaac4
#
#########################################

#!/bin/bash

# include config files
source ../configs/parameters.cfg

main()
{
	# so Isaac4 doesn't need to build an index first.
	#
	# you can build a metadata file first from the reference
	# which apparently reduces processing time, but it is not
	# required.
	#
	# doesn't seem to be working at the moment...
	# ...I think its memory issue. Isaac may only be for high 
	# powered servers. OK. 

	species=lambda-virus
	threads=4

	#sort_reference

	paired_end_align
}


#
# functions
#
sort_reference()
{
	${isaac}/isaac-sort-reference \
		-g ${dat_dir}/references/${species}.fa \
		-o ${dat_dir}/indices
		
}

paired_end_align()
{
	# general usage:
	# $ isaac-align -r <reference> -b <base calls> -m <memory limit> [optional arguments]
	${isaac}/isaac-align \
		-r ${dat_dir}/references/${species}.fa \
		-b ${dat_dir}/reads/${species}-short \
		-m 60 \
		-j 6 \
		-f fastq \
		-t ${tmp_dir} \
		--verbosity=3 \
		-o ${dat_dir}/aligned
}

#
#  run main
#
main
