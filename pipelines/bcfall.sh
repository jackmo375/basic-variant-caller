#
#  bcfall
#	* NGS variant calling pipeline
#	* using bcftools and bwa aligner.
#	* See:
#	* http://www.htslib.org/workflow/
#	* for notes on best practices.
#
#	Kevin Esoh
#
#########################################

#!/usr/bin/env bash

source ../includes/locations.sh
source ../includes/utilities.sh
source ${pip_dir}/modules/checkparams.mod.sh
source ${pip_dir}/modules/geneMapNGS.mod.sh


workflow() {
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]})

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	custom_call fq "checking read file quality..."

	custom_call trim "trimming read files..."

	custom_call bmap "mapping reads to the reference with bwa..."

	custom_call bcfcall "calling variants with bcftools..."

}

#
#  run workflow
#
workflow "$@"