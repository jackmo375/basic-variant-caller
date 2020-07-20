#
#  geneMapNGS
#	*NGS variant calling pipeline*
#
#	Kevin Esoh
#
#######################

#!/usr/bin/env bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() { argv=("$@")

	ref=NULL
	adap=NULL
	dname="$(readlink -f `pwd`)/"
	leadx=0
	trailx=0
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

	species=$(value_from_json $input_json '.species')
	t=$(value_from_json $input_json '.threads')
	read_type=$(value_from_json $input_json '.read_type')

	ref=${ref_dir}/${species}.fa
	dname=${rds_dir}/
	meta=${rds_dir}/${species}.txt

	##start

	##custom_call pfq "checking read file quality..."

	##custom_call pbmap "aligning reads with bwa..."

	custom_call pgmap "GATKv4 BWA Mapping/Alignment..."
	
}


#
#  tasks/functions
#
source ./geneMapNGS.tasks.sh


#
#  run workflow
#
workflow "$@"
