#
#	SIMPLE:
#		testing a single pipeline with simulations
#
#	*compares truth and estimated vcf files
#	for simulated SNP+indel data based on an input
#	reference sequence, and output the Jaccard index
#	for the truth and estimated vcf files.*
#
###############################################

source ../configs/parameters.cfg

workflow() {  argv=("$@")

	test_suite=simple
	input_json=${argv[0]}
	output_label=$(basename $input_json .json)

	pipeline=$(value_from_json $input_json '.pipeline')
	pipelineArgs=$(value_from_json $input_json '.pipelineArgs')

	input_json_pipeline=${pip_dir}/$pipeline.$pipelineArgs.json
	species=$(value_from_json $input_json_pipeline '.species')
	read_type=$(value_from_json $input_json_pipeline '.read_type')

	# 1. simulate input data
	custom_call simulate_reads "simulating reads..."

	# 2. run pipeline and compute the runtime
	START=$(date +%s.%N)
	${pip_dir}/$pipeline.sh $input_json_pipeline
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)

	# 3. process pipeline outputs
	custom_call compare_truth_est_vcf "comparing the simulation's truth vcf file against the variant caller output..."

	custom_call jaccard_index "computing jaccard index..."

	echo "pipeline runtime: $DIFF seconds"

}


#
#  tasks
#
simulate_reads()
{
	perl ${simulate} \
		--ref=${ref_dir}/${species}.fa \
		--prefix=${rds_dir}/${species} \
		--input $input_json \
		--output_vcf ${vcf_dir}/$pipeline.$pipelineArgs.truth.vcf

	$bcftools sort \
		${vcf_dir}/$pipeline.$pipelineArgs.truth.vcf \
		-o ${vcf_dir}/$pipeline.$pipelineArgs.truth.sorted.vcf
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
		${vcf_dir}/$pipeline.$pipelineArgs.raw.g.vcf \
		-Oz \
		-o ${vcf_dir}/$pipeline.$pipelineArgs.raw.g.vcf.gz

	$bcftools view \
		${vcf_dir}/$pipeline.$pipelineArgs.truth.sorted.vcf \
		-Oz \
		-o ${vcf_dir}/$pipeline.$pipelineArgs.truth.sorted.vcf.gz

	echo "building indices for the zipped vcf files"
	$bcftools index ${vcf_dir}/$pipeline.$pipelineArgs.raw.g.vcf.gz
	$bcftools index ${vcf_dir}/$pipeline.$pipelineArgs.truth.sorted.vcf.gz

	echo "computing the set difference and intersection vcf files"
	$bcftools isec \
		${vcf_dir}/$pipeline.$pipelineArgs.raw.g.vcf.gz \
		${vcf_dir}/$pipeline.$pipelineArgs.truth.sorted.vcf.gz \
		-p ${vcf_dir}/${output_label}	
}

jaccard_index() {
	# compute the jaccard index from the output of
	# bcftools isec
	# assumes:
	#	0000, 0001 are set difference files, and
	#	0002, 0003 are intersection files

	prefix=${output_label}

	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0000.vcf > ${tmp_dir}/$prefix.txt
	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0001.vcf >> ${tmp_dir}/$prefix.txt
	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0002.vcf >> ${tmp_dir}/$prefix.txt
	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0003.vcf >> ${tmp_dir}/$prefix.txt

	python $jaccard \
		--input ${tmp_dir}/$prefix.txt
}


#
#  run workflow
#
workflow "$@"

