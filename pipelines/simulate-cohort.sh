#
#  SIMULATE COHORT - simulation pipeline
#
#	* generating a cohort of artificial reads 
#	* from an input reference sequence
#
#	Jack Morrice
#
###################################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${pip_dir}/modules/checkparams.mod.sh

workflow() { 
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]} ['tmp_prefix']=${argv[2]} )

	custom_call check_input_json "checking simulation input json file was provided..."

	custom_call initialize_inputs_hash "initializing simulation input parameter values..."

	custom_call simulate_cohort_reads "simulating reads for the cohort..."

}


#
#  tasks
#
initialize_inputs_hash() {
	local status=0

	# 1. set default parameter values
	printf '  setting default parameter values...'
	inputs["cohort_id"]="c_0"
	inputs["ref"]=NULL
	inputs["read_type"]="pe"; # pe: paried end, sr: single read
	inputs["n_samples"]=1
	inputs["threads"]=1
	echo '...done'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.cohort_id'	inputs["cohort_id"]
	value_from_json ${inputs["input_json"]} '.ref_base'		inputs["ref_base"] && inputs["ref"]=${ref_dir}/${inputs["ref_base"]}
	value_from_json ${inputs["input_json"]} '.read_type'	inputs["read_type"]
	value_from_json ${inputs["input_json"]} '.n_samples'	inputs["n_samples"]
	value_from_json ${inputs["input_json"]} '.threads'		inputs["threads"]
	echo '...done'

	# 3. check that inputs make sense
	printf '  checking that parameter values make sense...'
	check_id "cohort" ${inputs["cohort_id"]} || { echo 'invalid choice of cohort id'; status=1; }
	check_ref || status=1
	check_int ${inputs["n_samples"]} n_samples || status=1
	[[ $status == 0 ]] && echo '...done'

	# 4. set up logging information
	set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }

	return $status
}

simulate_cohort_reads() {

	_mutate_reference || return 1

	_generate_reads || return 1

	_sort_gvcf_file 'truth' || return 1
}

_mutate_reference() {

	printf '  mutating the reference...'
	$simulate MutateReference \
		--input_ref ${inputs["ref"]} \
		--output_ref ${inputs['tmp_prefix']}${inputs["cohort_id"]}.fa \
		--input_json ${inputs["input_json"]} \
		--output_vcf ${vcf_dir}/${inputs["cohort_id"]}.truth.g.vcf \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& echo "...done."
}

_generate_reads() {
	local \
		s \
		prefix \
		sample_log_file \
		option_string \
		log_file_string

	echo "# sample file for cohort: ${inputs["cohort_id"]}" > ${rds_dir}/${inputs["cohort_id"]}.txt
	for s in $(seq 1 1 ${inputs["n_samples"]}); do
		prefix="${inputs["cohort_id"]}.s_$s"

		# create a new log file for each sample if they don't already exist:
		sample_log_file=${inputs["log_prefix"]}${prefix}.log
		[[ -s $sample_log_file ]] || > $sample_log_file

		echo -e "${prefix}.raw_1P.fq\t${prefix}.raw_2P.fq\ts_$s\tILLUMINA" >> ${rds_dir}/${inputs["cohort_id"]}.txt
	done

	option_string="GenerateReads \
		--input_ref ${inputs['tmp_prefix']}${inputs["cohort_id"]}.fa \
		--reads_prefix "${rds_dir}/"${inputs["cohort_id"]}.{3}.raw \
		--input_json ${inputs["input_json"]}"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  simulating reads..."
	run_in_parallel \
		$simulate \
		${rds_dir}/${inputs["cohort_id"]}.txt \
		"${option_string}" \
		${inputs["threads"]} \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }
}

_sort_gvcf_file() {
	local input_gvcf_stage=$1

	printf '  sorting truth gvcf files...'
	$bcftools sort \
		${vcf_dir}/${inputs["cohort_id"]}.${input_gvcf_stage}.g.vcf \
		-o ${vcf_dir}/${inputs["cohort_id"]}.${input_gvcf_stage}.sorted.g.vcf \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }
}


#
#  run workflow
#
workflow "$@"