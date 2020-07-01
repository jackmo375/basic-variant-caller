#
#	testing the gatk variant caller pipeline
#
#	we are using this workflow for inspiration:
#	https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.wdl
#
###############################################

source ../configs/parameters.cfg

workflow() {

	species=lambda-virus
	threads=1

	printf "${green}xx${nc} initiating read alignment and variant calling"; echo

	# if the right index does not exist, build it:
	build_index

	# align reads to the reference:
	single_read_align

	# convert sam files to bam files for use in variant caller:
	sam_to_bam
}


#
#  tasks
#
build_index()
{
	# option -a specifies the algorithm, with arguments:
	# 
	# 1. 'is' is for short sequences (like virus genomes)
	# 2. 'bwtsw' is for long sequences (like human genomes)


	if [[ ! -f "${dat_dir}/indices/${species}.bwa.amb" ]]; then
		printf "${green}xxxx${nc} building bwa index..."; echo
		${bwa}/bwa index \
			-p ${dat_dir}/indices/${species}.bwa \
			-a is \
			${dat_dir}/references/${species}.fa
	else
		printf "${green}xxxx${nc} skipping index build as index already exists..."; echo
	fi
	
}

single_read_align()
{
	printf "${green}xxxx${nc} aligning single reads with bwa..."; echo

	# basic usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	# <idxbase> :: the base of the index file names
	#				i.e. the names *with* full path, but no extension
	# <in1.fq>  :: input file of sequence reads in fastq/fasta format
	# [in2.fq]  :: (optional) input file of mates for the first read file, 
	#				if the input data is paired-end 

	${bwa}/bwa mem \
		${dat_dir}/indices/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		> ${sam_dir}/${species}.bwa.sr.sam
}


sam_to_bam() {
	printf "${green}xxxx${nc} converting sam files to bam format..."; echo

	${samtools} view \
		-S \
		-b ${sam_dir}/${species}.bwa.sr.sam \
		> ${bam_dir}/${species}.bwa.sr.bam
}

#
#  execute workflow
#
workflow
