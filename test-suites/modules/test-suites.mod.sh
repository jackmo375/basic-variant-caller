#
#	TEST SUITES - bash module
#
#	* module of bash functions for use in 
#	* constructing generic pipeline test
#	* suites. 
#
#	Jack Morrice
#
###########################################

#  summary:
#	1. initialize_inputs_hash()
#	2. ...

initialize_inputs_hash() {
	local \
		cohort_id_pi \
		cohort_id_si \
		n_samples \
		status=0

	# 1. set default parameter values
	printf 'setting default parameter values...'
	inputs["runs_per_pipeline"]=1
	inputs["cohort_id"]=NULL
	inputs["simulate_id"]=simulate-cohort
	inputs["simulate_inputs_id"]=basic
	inputs["pipeline_id"]=gatkall
	inputs["pipeline_inputs_id"]=basic
	echo 'done'

	# 2. update parameters with arguments from the input json file
	printf 'updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.runs_per_pipeline'	inputs["runs_per_pipeline"]
	value_from_json ${inputs["input_json"]} '.cohort_id'			inputs["cohort_id"]
	value_from_json ${inputs["input_json"]} '.simulate_id'			inputs["simulate_id"]
	value_from_json ${inputs["input_json"]} '.simulate_inputs_id'	inputs["simulate_inputs_id"]
	value_from_json ${inputs["input_json"]} '.pipeline_id'			inputs["pipeline_id"]
	value_from_json ${inputs["input_json"]} '.pipeline_inputs_id'	inputs["pipeline_inputs_id"]
	echo 'done'

	# 3. check that inputs make sense
	printf 'checking that parameter values make sense...'
	check_int ${inputs["runs_per_pipeline"]} runs_per_pipeline || status=1
	check_id "cohort" ${inputs["cohort_id"]} || status=1
	## are all the input cohort ids the same?
	value_from_json \
		${pip_dir}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json \
		'.cohort_id' \
		cohort_id_si
	value_from_json ${pip_dir}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json \
		'.cohort_id' \
		cohort_id_pi
	[ "${inputs["cohort_id"]}" == "$cohort_id_si" -a "$cohort_id_si" == "$cohort_id_pi" ] \
		|| { echo 'test suite ERROR: cohort ids do not match, please check all input json files and try again'; status=1; }
	[ -s ${pip_dir}/${inputs["simulate_id"]}.sh ] || { echo "test suite ERROR: ${inputs["simulate_id"]}.sh file not found or is empty"; status=1; }
	[ -s ${pip_dir}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json ] || { echo "test suite ERROR: ${pip_dir}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json file not found or is empty"; status=1; }
	[ -s ${pip_dir}/${inputs["pipeline_id"]}.sh ] || { echo "test suite ERROR: ${inputs["pipeline_id"]}.sh file not found or is empty"; status=1; }
	[ -s ${pip_dir}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json ] || { echo "test suite ERROR: ${pip_dir}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json file not found or is empty"; status=1; }

	## still need to write checks for the other input parameters...
	[[ $status == 0 ]] && echo 'done'

	# 4. set up logging information
	set_up_log_directory || { echo 'seting up log directory failed'; status=1; }

	return $status
}

compare_truth_est_vcf()
{
	# compare truth and estimated vcf files
	# truth:     vcf file produced by the simulator
	# estimated: vcf file produced by the variant caller

	# need to zip the vcf files, because bcf tools only accepts
	# gzipped inputs

	echo "zipping vcf files"
	$bcftools view \
		${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf \
		-Oz \
		-o ${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf.gz

	$bcftools view \
		${vcf_dir}/${inputs["cohort_id"]}.truth.sorted.vcf \
		-Oz \
		-o ${vcf_dir}/${inputs["cohort_id"]}.truth.sorted.vcf.gz

	echo "building indices for the zipped vcf files"
	$bcftools index ${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf.gz
	$bcftools index ${vcf_dir}/${inputs["cohort_id"]}.truth.sorted.vcf.gz

	echo "computing the set difference and intersection vcf files"
	$bcftools isec \
		${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf.gz \
		${vcf_dir}/${inputs["cohort_id"]}.truth.sorted.vcf.gz \
		-p ${vcf_dir}/${inputs["cohort_id"]}.isec	
}

jaccard_index() {
	local \
		tmp_prefix=${tmp_dir}/$(random_id)_ \
		__resultvar=$1 \
		myresult='NULL'

	# compute the jaccard index from the output of
	# bcftools isec
	# assumes:
	#	0000, 0001 are set difference files, and
	#	0002, 0003 are intersection files

	$gatk CountVariants \
		-V ${vcf_dir}/${inputs["cohort_id"]}.isec/0000.vcf > ${tmp_prefix}${inputs["cohort_id"]}.isec
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs["cohort_id"]}.isec/0001.vcf >> ${tmp_prefix}${inputs["cohort_id"]}.isec
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs["cohort_id"]}.isec/0002.vcf >> ${tmp_prefix}${inputs["cohort_id"]}.isec
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs["cohort_id"]}.isec/0003.vcf >> ${tmp_prefix}${inputs["cohort_id"]}.isec

	myresult=$(python $jaccard \
		--input ${tmp_prefix}${inputs["cohort_id"]}.isec)

	eval $__resultvar="'$myresult'"

}