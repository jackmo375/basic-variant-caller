
source ../configs/parameters.cfg

workflow() {
	clean_dir "$bam_dir/*" "bam"
	clean_dir "$sam_dir/*" "sam"
	clean_dir "$vcf_dir/*" "vcf"
	clean_dir "$tmp_dir/*" "temp"
	clean_dir "$dat_dir/references/*.fai" "reference index"
	clean_dir "$dat_dir/references/*.dict" "reference dict"
	clean_dir "$dat_dir/indices/*" "index"
}

#
#  functions
#
clean_dir() {	
	rm -rf $1 && echo "$2 folder clean"
}

workflow
