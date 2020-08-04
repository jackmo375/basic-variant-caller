#
#	SIMPLE:
#		testing a single pipeline with simulations of cohorts
#
#	*compares truth and estimated vcf files
#	for simulated SNP+indel data based on an input
#	reference sequence, and output the Jaccard index
#	for the truth and estimated vcf files.*
#
################################################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${pip_dir}/modules/checkparams.mod.sh
source ${tst_dir}/modules/test-suites.mod.sh


workflow () {
	local \
		argv=("$@") \
		jacc_value='NULL'

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]})

	custom_call check_input_json "checking input json file was provided..."

	custom_call initialize_inputs_hash "initializing input parameter values..."
	
	custom_call jaccard_index "computing jaccard index..." \
		jacc_value

	echo "Jaccard value is: $jacc_value"
	
}


#
#  tasks
#


#
#  run workflow
#
workflow "$@"