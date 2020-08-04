#  clean-data
#
#	* clean all intermediate and final PIPELINE output data
#
#	Jack Morrice
#
###############################################################

#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { 
	local argv=("$@")
	
	echo 'cleaning fastqc files' && rm ${fqc_dir}/* 2> /dev/null
	echo 'cleaning trimmed read files' && rm ${rds_dir}/*.paired* ${rds_dir}/*.unpaired* 2> /dev/null
	echo 'cleaning sam files' && rm ${sam_dir}/* 2> /dev/null
	echo 'cleaning bam files' && rm ${bam_dir}/* 2> /dev/null
	echo 'cleaning bam files' && rm -r ${vcf_dir}/* 2> /dev/null
	echo 'cleaning temporary files' && rm ${tmp_dir}/* 2> /dev/null
	#echo 'cleaning log files' && rm -r ${log_dir}/* 2> /dev/null
}

#
#  tasks
#

#
#  run workflow
#
workflow "$@"

