#
#  simulate-cohort
#	* generating a cohort of artificial reads from a reference
#
#	Jack Morrice
#
#################################################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${pip_dir}/modules/checkparams.mod.sh

workflow() { 
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} )

	custom_call check_input_json "checking input json file was provided..."

	custom_call initialize_inputs_hash "initializing input parameter values..."

	custom_call simulate_cohort_reads "simulating reads for the cohort..."

}


#
#  tasks
#
initialize_inputs_hash() {
	local status=0

	# 1. set default parameter values
	printf 'setting default parameter values...'
	inputs["cohort_id"]="c_0"
	inputs["ref"]=NULL
	inputs["read_type"]="pe"; # pe: paried end, sr: single read
	inputs["n_samples"]=1
	echo 'done'

	# 2. update parameters with arguments from the input json file
	printf 'updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.cohort_id'	inputs["cohort_id"]
	value_from_json ${inputs["input_json"]} '.ref_base'		inputs["ref_base"] && inputs["ref"]=${ref_dir}/${inputs["ref_base"]}
	value_from_json ${inputs["input_json"]} '.read_type'	inputs["read_type"]
	value_from_json ${inputs["input_json"]} '.n_samples'	inputs["n_samples"]
	echo 'done'

	# 3. check that inputs make sense
	printf 'checking that parameter values make sense...'
	check_id "cohort" ${inputs["cohort_id"]} || { echo 'invalid choice of cohort id'; status=1; }
	check_ref || status=1
	check_int ${inputs["n_samples"]} n_samples || status=1
	[[ $status == 0 ]] && echo 'done'

	return $status
}

simulate_cohort_reads() {
	local \
		tmp_prefix=${tmp_dir}/$(random_id)_

	$simulate MutateReference \
		--input_ref ${inputs["ref"]} \
		--output_ref ${tmp_prefix}${inputs["cohort_id"]}.fa \
		--input_json ${inputs["input_json"]} \
		--output_vcf ${vcf_dir}/${inputs["cohort_id"]}.truth.vcf

	_generate_reads $tmp_prefix || return 1

	$bcftools sort \
		${vcf_dir}/${inputs["cohort_id"]}.truth.vcf \
		-o ${vcf_dir}/${inputs["cohort_id"]}.truth.sorted.vcf

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}

_generate_reads() {
	local \
		tmp_prefix=$1 \
		sample_ids \
		s \
		prefix \
		status=0

	sample_ids=$(seq 1 1 ${inputs["n_samples"]})
	echo "# sample file for cohort: ${inputs["cohort_id"]}" > ${rds_dir}/${inputs["cohort_id"]}.txt
	for s in ${sample_ids[@]}; do
		prefix="${inputs["cohort_id"]}.s_$s.reads"
		$simulate GenerateReads \
			--input_ref ${tmp_prefix}${inputs["cohort_id"]}.fa \
			--reads_prefix ${rds_dir}/$prefix \
			--input_json ${inputs["input_json"]} \
			|| { echo "ERROR: generating reads failed for sample id s_$s"; status=1; }
		echo -e "${prefix}_1.fq\t${prefix}_2.fq\ts_$s\tILLUMINA" >> ${rds_dir}/${inputs["cohort_id"]}.txt
	done

	return $status
}

#
#  run workflow
#
workflow "$@"