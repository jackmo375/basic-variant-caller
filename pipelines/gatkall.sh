#
#  gatkall
#	* NGS variant calling pipeline
#	* using the gatk 
#
#	Kevin Esoh
#
####################################

#!/usr/bin/env bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { local argv=("$@")
	
	declare -A inputs=( ["input_json"]=${argv[0]} )

	custom_call initialize_inputs_hash "initializing input parameter values..."

	#echo ${inputs["ref"]}; exit 0

	##custom_call pfq "checking read file quality..."

	##custom_call ptrim "trimming read files..."

	##custom_call gmap "performing Mapping/Alignment with GATKv4 and BWA ..."

	custom_call varcall "calling variants with gatk..."

}


#
#  tasks/functions
#
source ./geneMapNGS.tasks.sh

initialize_inputs_hash() {

	# 1. set default parameter values
	printf 'setting default parameter values...'
	inputs["meta"]=NULL
	inputs["threads"]=1
	inputs["trim_minlen"]=36
	inputs["trim_adap"]=NULL
	inputs["trim_leadx"]=0
	inputs["trim_trailx"]=0
	inputs["ref"]=NULL
	inputs["blist"]=NULL
	echo 'done'

	# 2. update parameters with arguments from the input json file
	printf 'updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.meta_base'   inputs["meta_base"] && inputs["meta"]=${rds_dir}/${inputs["meta_base"]}
	value_from_json ${inputs["input_json"]} '.threads'     inputs["threads"]
	value_from_json ${inputs["input_json"]} '.trim_minlen' inputs["trim_minlen"]
	value_from_json ${inputs["input_json"]} '.trim_leadx'  inputs["trim_leadx"]
	value_from_json ${inputs["input_json"]} '.trim_trailx' inputs["trim_trailx"]
	value_from_json ${inputs["input_json"]} '.ref_base'	   inputs["ref_base"] && inputs["ref"]=${ref_dir}/${inputs["ref_base"]}
	value_from_json ${inputs["input_json"]} '.blist' inputs["blist"]
	echo 'done'

	# 3. check that inputs make sense
	printf 'checking that parameter values make sense...'
	local status=0
	check_sample || status=1
	check_int ${inputs["threads"]} threads || status=1
	check_int ${inputs["trim_minlen"]} trim_minlen || status=1
	check_int ${inputs["trim_leadx"]} trim_leadx || status=1
	check_int ${inputs["trim_trailx"]} trim_trailx || status=1
	check_ref || status=1
	check_bwa_idx || status=1
	check_gatk_dict || status=1
	check_samtools_fai || status=1
	[[ $status == 0 ]] && echo 'done'

	return $status
}

#
#  run workflow
#
workflow "$@"