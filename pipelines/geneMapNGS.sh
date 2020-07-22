#
#  geneMapNGS
#	*NGS variant calling pipeline*
#
#	Kevin Esoh
#
####################################

#!/usr/bin/env bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { argv=("$@")

	# set input parameter defaults

	ref=NULL
	adap=NULL
	dname=${dat_dir}
	leadx=0
	trailx=0
	minlen=36
	t=1
	out="ngs"
	meta=NULL
	par=".par.txt"
	ks=NULL
	blist=NULL
	glist=NULL
	dbsnp=NULL
	ped=NULL

	input_json=${argv[0]}

	# update input parameters from json file

	ref=${ref_dir}/$(value_from_json $input_json '.ref')
	t=$(value_from_json $input_json '.threads')
	meta=${rds_dir}/$(value_from_json $input_json '.meta')
	minlen=$(value_from_json $input_json '.minlen')

	##start

	##custom_call pfq "checking read file quality..."

	##custom_call ptrim "trimming read files..."

	##custom_call pbmap "aligning reads with bwa..."

	custom_call gmap "GATKv4 BWA Mapping/Alignment..."
	
}


#
#  tasks/functions
#
source ./geneMapNGS.tasks.sh


#
#  run workflow
#
workflow "$@"
