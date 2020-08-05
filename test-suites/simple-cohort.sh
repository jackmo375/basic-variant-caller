#
#	SIMPLE COHORT - bash pipeline
#
#	* compares truth and estimated vcf files
#	* for simulated SNP+indel data based on an input
#	* reference sequence, and output the Jaccard index
#	* for the truth and estimated vcf files.
#
#	Jack Morrice
#
################################################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${tst_dir}/modules/test-suites.mod.sh
source ${pip_dir}/modules/checkparams.mod.sh

workflow () {
	local \
		argv=("$@") \
		jacc_value='NULL'

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]})

	custom_call check_input_json "checking test suite input json file was provided..."

	custom_call initialize_inputs_hash "initializing test suite input parameter values..."

	# 1. simulate reads:
	${pip_dir}/${inputs["simulate_id"]}.sh \
		${pip_dir}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json \
		${inputs["log_prefix"]} \
		|| { echo 'test suite ERROR: simulating step failed'; exit 1; }

	# 2. run the pipeline and time it:
	START=$(date +%s.%N)
	${pip_dir}/${inputs["pipeline_id"]}.sh \
		${pip_dir}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json \
		${inputs["log_prefix"]} \
		|| { echo 'test suite ERROR: pipeline step failed'; exit 1; }
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)

	# 3. compare truth and estimated gvcf files

	custom_call compare_truth_est_vcf "comparing the simulation's truth vcf file against the variant caller output..."
	
	custom_call jaccard_index "computing jaccard index..." \
		jacc_value
	
	echo "  Jaccard value: $jacc_value"
	echo "  pipeline runtime: $DIFF seconds"

}


#
#  tasks
#


#
#  run workflow
#
workflow "$@"