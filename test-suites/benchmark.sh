#
#  BENCHMARK - test suite
#
#	* comparing multiple pipelines
#	* via comparison plots of speed,
#	* memory consumption, and accuracy
#
#	Jack Morrice
#
#########################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${pip_dir}/modules/checkparams.mod.sh
source ${tst_dir}/modules/test-suites.mod.sh

workflow() {
	local 
		argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]})

	custom_call check_input_json "checking input json file was provided..."

	custom_call initialize_inputs_hash "initializing input parameter values..."

	custom_call run_benchmark_test "running pipeline benchmark test..."

	custom_call plot_benchmark_test "plotting benchmark test results..."

}


#
#  tasks
#
run_benchmark_test() {
	local \
		jacc_value='NULL'

	> ${dat_dir}/benchmark-runs.dat
	echo -e "run\tpipeline\tjaccard\truntime" >> ${dat_dir}/benchmark-runs.dat

	for i in $(seq 1 ${inputs["runs_per_pipeline"]}); do

		# clean data directories
		${tst_dir}/clean-data.sh

		# 1. simulate reads:
		${pip_dir}/${inputs["simulate_id"]}.sh \
			${pip_dir}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json \
			${inputs["log_prefix"]} \
			|| { echo 'test suite ERROR: simulating step failed'; exit 1; }

		for j in "gatkall" "bcfall"; do

			inputs["pipeline_id"]=$j
			echo "pipeline is: ${inputs["pipeline_id"]}"

			# clean data directories
			${pip_dir}/clean-data.sh

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
			
			echo "Jaccard value: $jacc_value"
			echo "pipeline runtime: $DIFF seconds"

			echo -e "$i\t${inputs["pipeline_id"]}\t${jacc_value}\t$DIFF" >> ${dat_dir}/benchmark-runs.dat

		done

	done
}

plot_benchmark_test() {

	Rscript $benchmarkR \
		${dat_dir}/benchmark-runs.dat \
		${tls_dir}/jaccard \
		${med_dir}/benchmark

}

#
#  run workfow
#
workflow "$@"
