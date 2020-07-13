#
#	comparing truth and estimated vcf files
#	for simulated SNP+indel data
#
#	we are using this workflow for inspiration:
#	https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.wdl
#
###############################################

source ../configs/parameters.cfg

workflow() {

	species=lambda_virus
	threads=8
	read_type=pe # pe: paried end; sr: single_read

	custom_call simulate_reads "simulating reads..."

	custom_call build_index "building alignment index..."

	custom_call build_varcall_index "building variant caller index..."

	custom_call build_varcall_dict "building variant caller dictionary"

	if [[ ${read_type} == 'pe' ]]; then 
		custom_call paired_end_align "aligning paired-end reads with bwa..."
	else
		custom_call single_read_align "aligning reads to the reference..."
	fi

	custom_call sam_to_bam "converting sam files to bam format..."

	custom_call sort_sam "sorting bam files..."

	custom_call mark_duplicates "marking duplicates..."

	custom_call validate_bam "validating bam files..."

	# for now we will skip the recalibration of base quality scores. 

	custom_call index_bam "indexing all input bam files..."

	custom_call call_variants "calling variants..."

	custom_call compare_truth_est_vcf "comparing the simulation's truth vcf file against the variant caller output..."

	custom_call jaccard_index "computing jaccard index..."
}


#
#  tasks
#
simulate_reads()
{
	perl ${simulate} \
		--ref=${ref_dir}/${species}.fa \
		--prefix=${rds_dir}/${species} \
		--input ./simulate-ini.json \
		--output_vcf ${vcf_dir}/${species}.truth.vcf

	$bcftools sort \
		${vcf_dir}/${species}.truth.vcf \
		-o ${vcf_dir}/${species}.truth.sorted.vcf
}

build_index()
{
	# option -a specifies the algorithm, with arguments:
	# 
	# 1. 'is' is for short sequences (like virus genomes)
	# 2. 'bwtsw' is for long sequences (like human genomes)


	if [[ ! -f "${dat_dir}/indices/${species}.bwa.amb" ]]; then
		$bwa index \
			-p ${dat_dir}/indices/${species}.bwa \
			-a is \
			${dat_dir}/references/${species}.fa
	else
		printf "${green}skipping index build as index already exists${nc}"; echo
	fi
	
}

build_varcall_index() 
{
	# the gatk haplotype variant caller requires a dictionary and an index 
	# are both built from the reference. this function builds the index 
	# which I think is different from the alignment index?
	#
	# http://www.htslib.org/doc/samtools-faidx.html
	# https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

	if [[ ! -f "${dat_dir}/references/${species}.fai" ]]; then
		$samtools faidx \
			${dat_dir}/references/${species}.fa \
			> ${dat_dir}/references/${species}.fai
	else
		printf "${green}skipping varcall index build as this index already exists${nc}"; echo
	fi

}

build_varcall_dict()
{
	# the gatk haplotype variant caller requires a dictionary and an index 
	# are both built from the reference. this function builds the dictionary 
	# which I think is different from the alignment index?
	#
	# https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
	# https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-

	if [[ ! -f "${dat_dir}/references/${species}.dict" ]]; then
		java -jar $picard CreateSequenceDictionary \
			R=${dat_dir}/references/${species}.fa \
			O=${dat_dir}/references/${species}.dict
	else
		printf "${green}skipping varcall dictionary build as this dictionary already exists${nc}"; echo
	fi

}

single_read_align()
{
	# basic usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	# <idxbase> :: the base of the index file names
	#				i.e. the names *with* full path, but no extension
	# <in1.fq>  :: input file of sequence reads in fastq/fasta format
	# [in2.fq]  :: (optional) input file of mates for the first read file, 
	#				if the input data is paired-end 
	#
	# gatk requires read group information. I am not sure exactly what
	# this means, but there is information here:
	# https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups

	$bwa mem \
		${dat_dir}/indices/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		-R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:LIB-DAD-1\tSM:DAD\tPI:200" \
		-t $threads \
		> ${sam_dir}/${species}.bwa.${read_type}.sam
}

paired_end_align()
{
	echo 'aligning paired-end reads with bwa...'

	${bwa} mem \
		${dat_dir}/indices/${species}.bwa \
		${dat_dir}/reads/${species}_1.fq \
		${dat_dir}/reads/${species}_2.fq \
		-R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:LIB-DAD-1\tSM:DAD\tPI:200" \
		-t $threads \
		> ${sam_dir}/${species}.bwa.${read_type}.sam

}

sam_to_bam()
{
	${samtools} view \
		-S \
		-b ${sam_dir}/${species}.bwa.${read_type}.sam \
		> ${bam_dir}/${species}.bwa.${read_type}.bam
}

sort_sam()
{
	# the sam files need to be sorted before we can mark any duplicates. 
	# this can be done with either picard or samtools.
	#
	# even though the picard function is called 'SortSam' it seems to 
	# actually sort bam files. That's fine. 

	java -jar ${picard} SortSam \
		I=${bam_dir}/${species}.bwa.${read_type}.bam \
		O=${bam_dir}/${species}.bwa.${read_type}.sorted.bam \
		SORT_ORDER=coordinate
}

mark_duplicates()
{
	# this can be done using either samtools or picard.
	# http://www.htslib.org/algorithms/duplicate.html
	# https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

	java -jar ${picard} MarkDuplicates \
		I=${bam_dir}/${species}.bwa.${read_type}.sorted.bam \
		O=${bam_dir}/${species}.bwa.${read_type}.marked.bam \
		M=${bam_dir}/${species}.bwa.${read_type}.marked_dup_metrics.txt

}

validate_bam()
{
	java -jar $picard ValidateSamFile \
		I=${bam_dir}/${species}.bwa.${read_type}.marked.bam \
		MODE=SUMMARY
}

index_bam()
{
	$samtools index \
		${bam_dir}/${species}.bwa.${read_type}.marked.bam

}

call_variants()
{
	# I think we need to perform variant calling for a sample ensemble?
	# doesn't seem to work for a single sample? Not sure...

	$gatk --java-options "-Xmx4g" HaplotypeCaller  \
		--native-pair-hmm-threads $threads \
		-R ${dat_dir}/references/${species}.fa \
		-I ${bam_dir}/${species}.bwa.${read_type}.marked.bam \
		-O ${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf
}

compare_truth_est_vcf()
{
	# compare truth and estimated vcf files
	# truth:     vcf file produced by the simulator
	# estimated: vcf file produced by the variant caller

	# need tozip the vcf files, because bcf tools only accepts
	# gzipped inputs

	echo "zipping vcf files"
	$bcftools view \
		${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf \
		-Oz \
		-o ${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf.gz

	$bcftools view \
		${vcf_dir}/${species}.truth.sorted.vcf \
		-Oz \
		-o ${vcf_dir}/${species}.truth.sorted.vcf.gz

	echo "building indices for the zipped vcf files"
	$bcftools index ${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf.gz
	$bcftools index ${vcf_dir}/${species}.truth.sorted.vcf.gz

	echo "computing the set difference and intersection vcf files"
	$bcftools isec \
		${vcf_dir}/${species}.bwa.${read_type}.raw.g.vcf.gz \
		${vcf_dir}/${species}.truth.sorted.vcf.gz \
		-p ${vcf_dir}/compare	
}

jaccard_index() {
	# compute the jaccard index from the output of
	# bcftools isec
	# assumes:
	#	0000, 0001 are set difference files, and
	#	0002, 0003 are intersection files

	prefix=compare

	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0000.vcf > ${tmp_dir}/$prefix.txt
	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0001.vcf >> ${tmp_dir}/$prefix.txt
	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0002.vcf >> ${tmp_dir}/$prefix.txt
	$gatk CountVariants \
		-V ${vcf_dir}/$prefix/0003.vcf >> ${tmp_dir}/$prefix.txt

	python $jaccard \
		--input ${tmp_dir}/$prefix.txt
}


#
#  run workflow
#
workflow
