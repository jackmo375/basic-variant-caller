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
source ${pip_dir}/modules/checkparams.mod.sh
source ${pip_dir}/modules/geneMapNGS.mod.sh


workflow() {
	local \
		argv=("$@") \
	
	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]})

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	custom_call fq "checking read file quality..."

	custom_call trim "trimming read files..."

	custom_call gmap "performing Mapping/Alignment with GATKv4 and BWA ..."

	custom_call varcall "calling variants with gatk..."

}


#
#  run workflow
#
workflow "$@"