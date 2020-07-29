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

workflow() {
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} )

	custom_call initialize_inputs_hash "initializing input parameter values..."

	##custom_call fq "checking read file quality..."

	##custom_call trim "trimming read files..."

	custom_call bmap "map reads to the reference with bwa..."

	custom_call bcfcall "calling variants with bcftools..."

}

#
#  tasks/functions
#
source ./geneMapNGS.tasks.sh

#
#  run workflow
#
workflow "$@"