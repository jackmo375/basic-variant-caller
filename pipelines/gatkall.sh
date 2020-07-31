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

	custom_call initialize_inputs_hash "initializing input parameter values..."

	##custom_call fq "checking read file quality..."

	custom_call ptrim "trimming read files..."
	exit 0

	custom_call gmap "performing Mapping/Alignment with GATKv4 and BWA ..."

	custom_call varcall "calling variants with gatk..."

}


#
#  run workflow
#
workflow "$@"