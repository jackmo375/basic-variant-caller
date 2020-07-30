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

workflow () {
	local argv=("$@")

	testsuite_id=ts_1

	cohort_id=c_1

	simulate_id=simulate-cohort
	simulate_inputs_id=basic
	pipeline_id=gatkall
	pipeline_inputs_id=basic

	# 1. simulate reads:
	${pip_dir}/${simulate_id}.sh \
		${pip_dir}/${simulate_id}.${simulate_inputs_id}.json \
		|| { echo 'test suite ERROR: simulating step failed'; exit 1; }

	# 2. run the pipeline and time it:
	START=$(date +%s.%N)
	${pip_dir}/${pipeline_id}.sh \
		${pip_dir}/${pipeline_id}.${pipeline_inputs_id}.json \
		|| { echo 'test suite ERROR: pipeline step failed'; exit 1; }
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)

	# 3. compare truth and estimated gvcf files
	custom_call compare_truth_est_vcf "comparing the simulation's truth vcf file against the variant caller output..."
	custom_call jaccard_index "computing jaccard index..."
	echo "pipeline runtime: $DIFF seconds"

}


#
#  tasks
#
compare_truth_est_vcf()
{
	# compare truth and estimated vcf files
	# truth:     vcf file produced by the simulator
	# estimated: vcf file produced by the variant caller

	# need to zip the vcf files, because bcf tools only accepts
	# gzipped inputs

	echo "zipping vcf files"
	$bcftools view \
		${vcf_dir}/${cohort_id}.genotyped.g.vcf \
		-Oz \
		-o ${vcf_dir}/${cohort_id}.genotyped.g.vcf.gz

	$bcftools view \
		${vcf_dir}/${cohort_id}.truth.sorted.vcf \
		-Oz \
		-o ${vcf_dir}/${cohort_id}.truth.sorted.vcf.gz

	echo "building indices for the zipped vcf files"
	$bcftools index ${vcf_dir}/${cohort_id}.genotyped.g.vcf.gz
	$bcftools index ${vcf_dir}/${cohort_id}.truth.sorted.vcf.gz

	echo "computing the set difference and intersection vcf files"
	$bcftools isec \
		${vcf_dir}/${cohort_id}.genotyped.g.vcf.gz \
		${vcf_dir}/${cohort_id}.truth.sorted.vcf.gz \
		-p ${vcf_dir}/${cohort_id}.isec	
}

jaccard_index() {
	# compute the jaccard index from the output of
	# bcftools isec
	# assumes:
	#	0000, 0001 are set difference files, and
	#	0002, 0003 are intersection files

	prefix=${cohort_id}.isec

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