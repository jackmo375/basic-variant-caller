
source ../configs/parameters.cfg

workflow() {

	species=lambda_virus
	threads=1
	read_type=pe # pe: paried end; sr: single_read

	#custom_call compare_truth_est_vcf "comparing the simulation's truth vcf file against the variant caller output..."

	custom_call jaccard_index "computing jaccard index..."
}

compare_truth_est_vcf()
{
	# compare truth and estimated vcf files
	# truth:     vcf file produced by the simulator
	# estimated: vcf file produced by the variant caller

	# need tozip the vcf files, because bcf tools only accepts
	# gzipped inputs

	outputDir=${vcf_dir}/compare

	echo "zipping vcf files"
	$bcftools view \
		${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf \
		-Oz \
		-o ${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf.gz \
	$bcftools view \
		${vcf_dir}/${species}.truth.sorted.vcf \
		-Oz \
		-o ${vcf_dir}/${species}.truth.sorted.vcf.gz \
	
	echo "building indices for the zipped vcf files"
	$bcftools index ${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf.gz
	$bcftools index ${vcf_dir}/${species}.truth.sorted.vcf.gz

	echo "computing the set difference and intersection vcf files"
	$bcftools isec \
		${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf.gz \
		${vcf_dir}/${species}.truth.sorted.vcf.gz \
		-p $outputDir
}

jaccard_index() {
	# compute the jaccard index from the output of
	# bcftools isec
	# assumes:
	#	0000, 0001 are set difference files, and
	#	0002, 0003 are intersection files

	prefix=compare

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
workflow
