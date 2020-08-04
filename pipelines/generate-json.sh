#
#	GENERATE JSON
#
#	* generate an input json file series
#	* to input to pipelines (for testing
#	* the variation of pipeline output
#	* with input parameters). 
#
#	Jack Morrice
#
##########################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${pip_dir}/modules/checkparams.mod.sh

workflow() {
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} )

	custom_call check_input_json "checking input json file is valid"

	python ${generatejson} \
		--input_json ${inputs["input_json"]} \
		--key 'cohort_id' \
		--new_value 'c_2' \
		--output_json 'out.json'
}

#
#  tasks
#
initialize_inputs_hash() {
	echo 'test'
}

#
#  run workflow
#
workflow "$@"