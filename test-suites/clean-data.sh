#  clean-data
#
#	* clean all intermediate and final TEST SUITE output data
#
#	Jack Morrice
#
###############################################################

#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { 
	local argv=("$@")

	${pip_dir}/clean-data.sh

	echo '  cleaning reads files' && rm ${rds_dir}/* 2> /dev/null
	echo '  cleaning truth vcf files' && rm ${vcf_dir}/*.truth.* 2> /dev/null
	echo '  cleaning isec vcf files' && rm -r ${vcf_dir}/*.isec* 2> /dev/null

}

#
#  tasks
#

#
#  run workflow
#
workflow "$@"