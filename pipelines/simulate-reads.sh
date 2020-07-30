#
#  simulate
#	*generating artificial reads from a reference sequence*
#	
#	Jack Morrice
#
##############################################################

source ../includes/locations.sh
source ../includes/utilities.sh


workflow() { 
	local argv=("$@")

	input_json=${argv[0]}
	output_label=$(basename $input_json .json)

	species=$(value_from_json $input_json '.species')
	read_type=$(value_from_json $input_json '.read_type')

	# simulate input data
	custom_call simulate_reads "simulating reads..."
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
		--output_vcf ${vcf_dir}/${output_label}.truth.vcf

	$bcftools sort \
		${vcf_dir}/${output_label}.truth.vcf \
		-o ${vcf_dir}/${output_label}.truth.sorted.vcf
}

#
#  run workflow
#
workflow "$@"
