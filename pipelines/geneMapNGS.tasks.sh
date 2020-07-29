#
#  function library for geneMapNGS pipeline
#
#	Kevin Esoh
#
#	contents:
#		1. command functions
#			*these functions are called at the top level
#			by the pipeline scripts*
#		2. check functions
#			*these functions check that the necessary 
#			files exist for the command functions to run
#			properly (called by the command functions 
#			themselves; not called at pipeline level).*
#		3. prep functions
#			*these functions prepare the user's project
#			directory for the command functions to run
#			(again, called by the command functions 
#			themselves; not called at pipeline level)*
#
#############################################################

function start() {
    echo -e """
               ===================================================================
               \e[38;5;43mGeneMAP NGS Pipeline			          GeneMAP (c) 2020\e[0m
               -------------------------------------------------------------------
               Argument:            Parameter
               --------             --------
               Fastq/SAM/BAM path:  $dname
               reference:           $ref
	       	   dbsnp:		        $dbsnp
               leading:             $leadx
               trailing:            $trailx
               PED file:            $ped
               threads:             $t
               adapter:             $adap
               sample file:         $meta
               BAM list:            $blist
               GVCF list:           $glist
               known sites:         $ks
               outFile:             ${out}.vcf.gz
               ===================================================================
               Starting NGS Pipeline. Please wait...
    """
}

initialize_inputs_hash() {
	local status=0

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
	inputs["glist"]=NULL
	inputs["ped"]=NULL
	inputs["dbsnp"]=NULL
	echo 'done'

	# 2. update parameters with arguments from the input json file
	printf 'updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.meta_base'   inputs["meta_base"] && inputs["meta"]=${rds_dir}/${inputs["meta_base"]}
	value_from_json ${inputs["input_json"]} '.threads'     inputs["threads"]
	value_from_json ${inputs["input_json"]} '.trim_minlen' inputs["trim_minlen"]
	value_from_json ${inputs["input_json"]} '.trim_leadx'  inputs["trim_leadx"]
	value_from_json ${inputs["input_json"]} '.trim_trailx' inputs["trim_trailx"]
	value_from_json ${inputs["input_json"]} '.ref_base'	   inputs["ref_base"] && inputs["ref"]=${ref_dir}/${inputs["ref_base"]}
	value_from_json ${inputs["input_json"]} '.blist' 	   inputs["blist"]
	value_from_json ${inputs["input_json"]} '.glist' 	   inputs["glist"]
	value_from_json ${inputs["input_json"]} '.ped' 	   	   inputs["ped"]
	value_from_json ${inputs["input_json"]} '.dbsnp' 	   inputs["dbsnp"]
	echo 'done'

	# 3. check that inputs make sense
	printf 'checking that parameter values make sense...'
	check_sample || status=1
	check_int ${inputs["threads"]} threads || status=1
	check_int ${inputs["trim_minlen"]} trim_minlen || status=1
	check_int ${inputs["trim_leadx"]} trim_leadx || status=1
	check_int ${inputs["trim_trailx"]} trim_trailx || status=1
	check_ref || status=1
	check_bwa_idx || status=1
	check_gatk_dict || status=1
	check_samtools_fai || status=1
	## still need to write checks for the other input parameters...
	[[ $status == 0 ]] && echo 'done'

	return $status
}

#
#  1. command functions
#
#  summary:
#	+ fq(), pfq()
#	+ trim(), ptrim()
#	+ bmap(), pbmap()
#	+ gmap(), pgmap()
#	+ indelreal(), pindelreal()
#	+ bqsr(), pbqsr()
#	+ emit_gvcfs(), pemit_gvcfs()
#	+ combinegvcfs()
#	+ genogvcfs()
#	+ varcall()
#	+ vqsr()
#	+ bcfcall()
#	

#  fq(), pfq()
#
#	*calls FastQC on read files in series, parallel*
#
#	input:
#		+ $meta
#	tools required:
#		+ $fastqc
#		+ gnu parallel (for pfq() only)
#	outputs:
#		+ fastq quality report files in ${fqc_dir}
#		+ forward_reverse.txt in ${tmp_dir} 
#
function fq() {

	local tmp_prefix=${tmp_dir}/$(random_id)_
	local id=${tmp_prefix}fastq.input.txt
	local odr="${fqc_dir}/"

	prep_fq $tmp_prefix || return 1

    while read -r line; do
        echo -e "running FastQC"
        $fastqc -t ${inputs["threads"]} $line -o $odr
    done < $id

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

}

function pfq() {
	local tmp_prefix=${tmp_dir}/$(random_id)_
	local id=${tmp_prefix}fastq.input.txt
	local odr="${fqc_dir}/"
	local n=$((50/${inputs["threads"]}))

	prep_fq $tmp_prefix || return 1

	echo -e "running FastQC"
	echo -e "Your jobs will be split across $n parallel threads"
	cat $id | \
		parallel --col-sep ' ' echo -e "-t ${inputs["threads"]} {1} $(if [ -n {2} ]; then echo {2}; fi) -o $odr" | \
		xargs -I input -P$n sh -c "$fastqc input"

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}

#  tim(), ptrim()
#
#	*clips unwanted reads in fastq files with trimmomatic*
#	inputs:
#		+ $meta
#		+ forward_reverse.txt in ${tmp_dir} 
#	tools required:
#		+ $trimmomatic
#	outputs:
#		+ trimmed (and zipped) read (fastq) files
#
function trim() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	prep_trim $tmp_prefix || return 1

	id=${tmp_prefix}trim.input.txt
	while read -r line; do
		echo -e "running Trimmomatic"
		java -jar $trimmomatic PE \
              -phred33 $line \
              $(if [[ ${inputs["trim_adap"]} != NULL ]]; then check_adapter; fi) \
              LEADING:${inputs["trim_leadx"]} \
              TRAILING:${inputs["trim_trailx"]} \
              SLIDINGWINDOW:4:15 \
              MINLEN:${inputs["trim_minlen"]} \
              -threads ${inputs["threads"]} \
              || { echo 'trimmomatic failed'; return 1; }
	done < $id
	
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}

function ptrim() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	prep_trim $tmp_prefix || return 1

	id=${tmp_prefix}trim.input.txt
	n=$((50/${inputs["threads"]}))
	echo -e "running Trimmomatic"
	echo -e "Your jobs will be split across $n parallel threads"
	cat $id | \
		parallel --col-sep ' ' echo "PE -phred33 {} $(if [[ ${inputs["trim_adap"]} != NULL ]]; then check_adapter; fi) LEADING:${inputs["trim_leadx"]} TRAILING:${inputs["trim_trailx"]} SLIDINGWINDOW:4:15 MINLEN:${inputs["trim_minlen"]} -threads ${inputs["threads"]}" | \
		xargs -I input -P$n sh -c "java -jar $trimmomatic input" \
		|| { echo 'trimmomatic failed'; return 1; }
	
	# remove temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}

#  bmap(), pbmap()
#
#	* Mapping/Alignment (BWA)
#	* this function does three things in sequence:
#	*	1. aligns all reads to the reference sequence with bwa mem
#	*		 to create sam files...
#	*	2. ...converts the sam files to bam files with $samtools view...
#	*	3. ...and then sorts the bam files into 'mapped' bam files,
#	*		and removes the unsorted (unmapped bam files).
#
#	inputs:
#		+ $t
#		+ $ref
#
#	tools required:
#		+ $bwa
#		+ $samtools
#
#	outputs:
#		+ aligned and sorted (mapped) sequences for each sample, 
#			saved as compressed bam files
#
function bmap() {

	echo -e "running BWA-BCFTOOLS Alignment/Mapping: Serial\n"

	_aligning_reads_to_reference || return 1

	_converting_sam_to_bam || return 1

	_sorting_bam_files || return 1

	_mark_duplicates || return 1

	_validate_bam_files || return 1
}

function pbmap() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample || return 1
	check_ref || return 1
	prep_map $tmp_prefix || return 1

	id=${tmp_prefix}align.input.txt

	# create file names
	awk '{print $4,"-O BAM -o",$4}' ${tmp_prefix}align.input.txt | \
	   sed 's/\.sam/\.bam/2' > ${tmp_prefix}sam2bam.input.txt
	awk '{print $5,$5}' ${tmp_prefix}sam2bam.input.txt | \
	   sed 's/\.bam/\.mapped.bam/1' > ${tmp_prefix}sortbam.input.txt
	n=$((50/$t))
	echo -e "running BWA-BCFTOOLS Alignment/Mapping: parallel\n"

	# align each sample's set of reads to the reference sequence
	cat ${tmp_prefix}align.input.txt | \
	   parallel --col-sep ' ' echo "mem -t $t $ref {}" | \
	   xargs -I input -P$n sh -c "$bwa input"

	# compress each sample's set of aligned reads to bam format
	cat ${tmp_prefix}sam2bam.input.txt | \
	   parallel --col-sep ' ' echo "view -h {}" | \
	   xargs -I input -P$n sh -c "$samtools input"

	# sort each sample's set of aligned reads 
	cat ${tmp_prefix}sortbam.input.txt | \
	   parallel --col-sep ' ' echo "sort -O BAM --reference $ref -@ $t -o {}" | \
	   xargs -I input -P$n sh -c "$samtools input"

	# then remove all unsorted bam files
	for sam in $(awk '{print $4}' ${tmp_prefix}align.input.txt); do
	   rm ${sam};
	done
	for bam in $(awk '{print $2}' ${tmp_prefix}sortbam.input.txt); do
	   rm ${bam};
	done

	# not sure...
	for i in out.vcf.gz aligned/*.sam; do
	   if [[ -e "${i}" ]]; then
		  rm $i;
	   fi;
	done

	# remove all temporary files from this call
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

}

#  gmap(), pgmap()
#
#	* GATKv4 BWA Mapping/Alignment
#	* its odd: only the parallel version marks duplicates
#	* and indexes the final bam files. Is this a mistake?
# 
#	inputs:
#		+ $meta
#		+ $t
#		+ $ref
#		+ a bwa index for $ref
#		+ a gatk dictionary for $ref
#
#	tools required:
#
#	outputs:
#		+ [gmap()]  ${bam_dir}/<sample_id>.merged.bam (mapped bam files for each sample)
#
function gmap() {

	echo "performing alignment; your jobs will run in serial"

	_converting_fastq_to_sam || return 1

	_aligning_reads_to_reference || return 1

	_converting_sam_to_bam || return 1

	_sorting_bam_files || return 1

	_merge_bam_alignments || return 1

	_mark_duplicates || return 1

	_validate_bam_files || return 1

}

_converting_fastq_to_sam() {
	local \
		sample_id \
		line \
		status=0

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$gatk FastqToSam \
			-F1 ${rds_dir}/$(echo $line | awk '{print $1}') \
			-F2 ${rds_dir}/$(echo $line | awk '{print $2}') \
			-SM $sample_id \
			-PL $(echo $line | awk '{print $4}') \
			-RG $sample_id \
			-O ${bam_dir}/${sample_id}.unmapped.bam \
			|| { echo 'FastqToSam step failed'; status=1; }

	done < ${inputs["meta"]}

	return $status
}

_aligning_reads_to_reference() {
	local \
		sample_id \
		line \
		status=0

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$bwa mem \
			-R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" \
			-t ${inputs["threads"]} \
			${inputs["ref"]} $(echo $line | awk -v d="${rds_dir}/" '{print d$1,d$2}') \
			-o ${sam_dir}/${sample_id}.sam \
			|| { echo 'bwa mem step failed'; status=1; }
	done < ${inputs["meta"]}

	return $status
}

_converting_sam_to_bam() {
	local \
		sample_id \
		line \
		status=0

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$samtools view \
			-h ${sam_dir}/${sample_id}.sam \
			-O BAM \
			-o ${bam_dir}/${sample_id}.bam \
			|| { echo 'sam to bam step failed'; status=1; }
	done < ${inputs["meta"]}

	return $status
}

_sorting_bam_files() {
	local \
		sample_id \
		line \
		status=0

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$samtools sort \
			-O BAM \
			--reference ${inputs["ref"]} \
			-@ ${inputs["threads"]} \
			-o ${bam_dir}/${sample_id}.sorted.bam \
			${bam_dir}/${sample_id}.bam \
			|| { echo 'sorting the bam files step failed'; status=1; }
	done < ${inputs["meta"]}

	return $status
}

_merge_bam_alignments() {
	local \
		sample_id \
		line \
		status=0

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$gatk MergeBamAlignment \
			-R ${inputs["ref"]} \
			-UNMAPPED ${bam_dir}/${sample_id}.unmapped.bam \
			-ALIGNED ${bam_dir}/${sample_id}.sorted.bam \
			-O ${bam_dir}/${sample_id}.merged.bam \
			|| { echo 'merge bam alignment step failed'; status=1; }
	done < ${inputs["meta"]}

	return $status
}

_mark_duplicates() {
	local \
		input_bam_stage='merged' \
		sample_id \
		line \
		status=0

	# check if merge_bam_alignment step was skipped:
	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')
		if [[ ! -s ${bam_dir}/${sample_id}.merged.bam ]]; then
			input_bam_stage='sorted'
			break
		fi
	done < ${inputs["meta"]}

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$gatk MarkDuplicates \
			-I ${bam_dir}/${sample_id}.${input_bam_stage}.bam \
			-O ${bam_dir}/${sample_id}.marked.bam \
			-M ${bam_dir}/${sample_id}.marked_dup_metrics.txt \
			|| { echo 'marking duplicates step failed'; status=1; }
	done < ${inputs["meta"]}

	return $status
}

_validate_bam_files() {
	local \
		sample_id \
		line \
		status=0

	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		sample_id=$(echo $line | awk '{print $3}')

		$gatk ValidateSamFile \
			-I ${bam_dir}/${sample_id}.marked.bam \
			-R ${inputs["ref"]} \
			--TMP_DIR ${tmp_dir}/ \
			-M SUMMARY \
			|| { echo 'error validating bam files'; status=1; }
	done < ${inputs["meta"]}

	return $status
}


function pgmap() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample || return 1
	check_ref || return 1
	check_bwa_idx || return 1
	check_gatk_dict || return 1 
	check_samtools_fai || return 1

	id="${tmp_prefix}metadat.txt"
	awk '{print $1,$2,$3,$4}' ${meta} | sed '/^#/d' > $id

	n=$((50/$t))
	echo -e "GATK-BWA Alignment and Mark Duplicates.\nYour jobs will be split across $n parallel threads\n"
	#--- Make unmapped BAM files from raw FASTQ files
	cat ${id} | parallel --col-sep ' ' echo "FastqToSam -F1 ${rds_dir}/{1} -F2 ${rds_dir}/{2} -SM {3} -PL {4} -RG {3} -O ${bam_dir}/{3}.unmapped.bam" | xargs -I input -P$n sh -c "$gatk input" \
			|| { echo 'FastqToSam step failed'; return 1; }
	# align each sample's set of reads to the reference
	cat ${id} | parallel --col-sep ' ' echo "mem -t $t $ref ${rds_dir}/{1} ${rds_dir}/{2} -o ${sam_dir}/{3}.mapped.sam" | xargs -I input -P$n sh -c "$bwa input" \
			|| { echo 'bwa mem step failed'; return 1; }
	# convert aligned sequence sam files to bam
	cat ${id} | parallel --col-sep ' ' echo "view -O BAM -h ${sam_dir}/{3}.mapped.sam -o ${bam_dir}/{3}.unsorted.mapped.bam" | xargs -I input -P$n sh -c "$samtools input" \
			|| { echo 'sam to bam step failed'; return 1; }
	# sort the bam files
	cat ${id} | parallel --col-sep ' ' echo "sort -O BAM --reference $ref -@ $t -o ${bam_dir}/{3}.mapped.bam ${bam_dir}/{3}.unsorted.mapped.bam" | xargs -I input -P$n sh -c "$samtools input" \
			|| { echo 'sorting the bam files step failed'; return 1; }
	# merge bam alignment files
	cat ${id} | parallel --col-sep ' ' echo "MergeBamAlignment -O ${bam_dir}/{3}.bam -R ${ref} -UNMAPPED ${bam_dir}/{3}.unmapped.bam -ALIGNED ${bam_dir}/{3}.mapped.bam" | xargs -I input -P$n sh -c "$gatk input" \
			|| { echo 'merge bam alignment step failed'; return 1; }

	#--- Mark duplicates
	cat ${id} | parallel --col-sep ' ' echo MarkDuplicates -I ${bam_dir}/{3}.bam -O ${bam_dir}/{3}_mkdups.bam -M ${bam_dir}/{3}_marked_dup_metrics.txt --REMOVE_DUPLICATES false | xargs -I input -P$n sh -c "$gatk input" \
			|| { echo 'failed to mark duplicates'; return 1; }
	cat ${id} | parallel --col-sep ' ' echo index -b ${bam_dir}/{3}_mkdups.bam | xargs -I input -P$n sh -c "$samtools input" \
			|| { echo 'failed to index bam files'; return 1; }

	# remove all temporary files
	{
		cat ${id} | parallel --col-sep ' ' rm ${bam_dir}/{3}.bam && \
		cat ${id} | parallel --col-sep ' ' rm ${sam_dir}/{3}.mapped.sam && \
		cat ${id} | parallel --col-sep ' ' rm ${bam_dir}/{3}*mapped.bam
	} || { echo 'failed to remove intermediate bam and sam files'; return 1; }
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

}

#  indelreal(), pindelreal()
#
#	* Indel Realignment 
#	* (THis will not be run if GATK is used since HaplotypeCaller 
#	* essentially does local rearrangements).
#	* THIS FUNCTION DOES NOT YET WORK IN THE NEW CONTEXT
#
function indelreal() {
       check_ref; check_bamlist
       #--- Indel realignment (This requires GATKv3.x. Point to your installation of it in 'gatk3_esoh' above)
       ###  According to GATK Best Practices, this step is not necessary in the new pipeline, as HaplotypeCaller does a good job  ###
       #awk '{print $1,$2,$3,$4}' ${meta} > metadat.txt
       id="bam.list"
       n=$((50/$t))
       mkdir -p realigned
       while read -r line; do
            gatk -T RealignerTargetCreator -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "rtc.ks.txt" -a -s "rtc.ks.txt" ]; then cat rtc.ks.txt; fi; fi) -I aligned/${line} -o realigned/${line/.bam/.intervals}
            gatk -T IndelRealigner -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "ir.ks.txt" -a -s "ir.ks.txt" ]; then cat ir.ks.txt; fi; fi) -I aligned/${line} -targetIntervals realigned/${line/.bam/.intervals} -o realigned/${line/.bam/.realigned.bam}
            samtools index -b aligned/${line/.bam/.realigned.bam}
       done < ${id}
       if [ -e "ir.ks.txt" -o -e "rtc.ks.txt" ]; then rm ir.ks.txt || rm rtc.ks.txt; fi
}

function pindelreal() {
       check_ref; check_bamlist
       #--- Indel realignment (This requires GATKv3.x. Point to your installation of it in 'gatk3_esoh' above)
       ###  According to GATK Best Practices, this step is not necessary in the new pipeline, as HaplotypeCaller does a good job  ###
       #awk '{print $1,$2,$3,$4}' ${meta} > metadat.txt
       id="bam.list"
       n=$((50/$t))
       mkdir -p realigned
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo -T RealignerTargetCreator -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "rtc.ks.txt" -a -s "rtc.ks.txt" ]; then cat rtc.ks.txt; fi; fi) -I aligned/{1}.bam -o realigned/{1}.intervals | xargs -I input -P$n sh -c "$gatk3 input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo -T IndelRealigner -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "ir.ks.txt" -a -s "ir.ks.txt" ]; then cat ir.ks.txt; fi; fi) -I aligned/{1}_mkdups.bam -targetIntervals realigned/{1}.intervals -o realigned/{1}.realigned.bam | xargs -I input -P$n sh -c "$gatk3 input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo index -b aligned/{1}.realigned.bam | xargs -I input -P$n sh -c "samtools input"
       if [ -e "ir.ks.txt" -o -e "rtc.ks.txt" ]; then rm ir.ks.txt || rm rtc.ks.txt; fi
}

#  bqsr(), pbqsr()
#
#	* Base Quality Score Recallibration (BQSR)
#
#	inputs:
#		+ bam list (optional)
#		+ The input read data whose base quality scores need to be assessed.
#		+ A database of known polymorphic sites to skip over.
#	tools:
#		+ gatk 4.0 (BaseRecalibrator, Apply BQSR)
#	outputs:
#		+ A GATK Report file with many tables
function bqsr() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_ref || return 1
	check_gatk_dict || return 1
	check_bamlist $tmp_prefix || return 1

	id=${tmp_prefix}bam.list
	n=$((50/$t))

	while read -r line; do
		echo $line
		return 1
		$gatk BaseRecalibrator \
			-I ${dname}/${line/.bam/} \
			-R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "bqsr.ks.txt" -a -s "bqsr.ks.txt" ]; then rm rtc.ks.txt ir.ks.txt; cat bqsr.ks.txt; fi; fi) \
			-O bqsr/${line/.bam/_recal_data.table} \
			|| { echo "BaseRecalibrator step failed"; return 1; }
		$gatk ApplyBQSR \
			-R $ref 
			-I ${dname}/${line/.bam/} 
			--bqsr-recal-file bqsr/${line/.bam/_recal_data.table} \
			-O aligned/${line/.bam/.mapped.bam} \
			|| { echo "ApplyBQSR step failed"; return 1; }
	done < ${id}
	if [ -e "${id}" ]; then rm ${id}; fi

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

}

function pbqsr() {
       if [[ "$dname" == NULL ]]; then
          echo -e "\n\e[38;5;1mERROR\e[0m: -p,--path not provided! Please specify path to BAM files"; 1>&2;
          return 1;
       fi
       check_ref; check_gatk_dict; check_bamlist
       id="$blist"
       n=$((50/$t))
       mkdir -p bqsr; mkdir -p aligned
       #--- Base quality score recalibration (BQSR)
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo BaseRecalibrator -I ${dname}/{1/}.bam -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "bqsr.ks.txt" -a -s "bqsr.ks.txt" ]; then rm rtc.ks.txt ir.ks.txt; cat bqsr.ks.txt; fi; fi) -O bqsr/{1/}_recal_data.table | xargs -I input -P$n sh -c "gatk input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo ApplyBQSR -R $ref -I ${dname}/{1/}.bam --bqsr-recal-file bqsr/{1/}_recal_data.table -O aligned/{1/}.mapped.bam | xargs -I input -P$n sh -c "gatk input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' rm ${dname}/{1/}.bam
       if [ -e "${id}" ]; then rm ${id}; fi
}

#  varcall
#
#	* Variant Calling with GATK (Single Cohort Joint) in Serial
#	* This function takes analysis ready bam files, and:
#	* 1. indexes them
#	* 2. calls variants with the gatk haplotype caller
#	* 3. ...
#	* the bam files are the output of the bqsr step of the 
#	* pipeline. 
#	inputs:
#
#	tools required:
#
#	outputs:
#
function varcall() {
	local \
		tmp_prefix=${tmp_dir}/$(random_id)_ \
		input_bam_stage='marked' \
		input_gvcf_stage='raw'

	prep_bamlist $tmp_prefix $input_bam_stage || return 1

	_index_bam_files ${tmp_prefix}bam.list || return 1

	_call_sample_variants ${tmp_prefix}bam.list || return 1

	prep_gvcflist $tmp_prefix $input_gvcf_stage || return 1

	_combine_sample_gvcfs ${tmp_prefix}gvcf.list || return 1

	_genotype_combined_gvcf || return 1

	# remove all temporary files from this call
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

}

_index_bam_files() {
	local \
		bamlist=$1 \
		line \
		sample_id \
		bam_base \
		bam_file \
		status=0

	while read -r line; do
		sample_id=$(echo $line | awk '{print $1}')
		bam_base=$(echo $line | awk '{print $2}')
		bam_file=${bam_dir}/$bam_base

		# 1. index the analysis-ready bam file
		printf "Creating SAMTOOLS index: $bam_base..."
		$samtools index \
			-b \
			-@ ${inputs["threads"]} \
			$bam_file \
			&& echo 'done' || { echo 'indexing $bam_base failed'; status=1; }

	done < $bamlist

	return $status
}

_call_sample_variants() {
	local \
		bamlist=$1 \
		sample_id \
		bam_base \
		bam_file \
		status=0

	while read -r line; do
		sample_id=$(echo $line | awk '{print $1}')
		bam_base=$(echo $line | awk '{print $2}')
		bam_file=${bam_dir}/$bam_base

		# 2. call variants
		$gatk HaplotypeCaller \
			-R ${inputs["ref"]} \
			-I $bam_file \
			-O ${vcf_dir}/${sample_id}.raw.g.vcf \
			$(if [[ ${inputs["ped"]} != NULL ]]; then echo -ped ${inputs["ped"]}; fi) \
			$(if [[ ${inputs["dbsnp"]} != NULL ]]; then echo --dbsnp ${inputs["dbsnp"]}; fi) \
			--lenient true \
			-ERC GVCF \
			|| { echo "HaplotypeCaller failed for $sample_id"; status=1; }

	done < ${tmp_prefix}bam.list

	return $status
}

_combine_sample_gvcfs() {
	local \
		gvcf_list=$1 \
		status=0

	$gatk CombineGVCFs \
		-R ${inputs["ref"]} \
		--arguments_file $gvcf_list \
		$(if [[ ${inputs["dbsnp"]} != NULL ]]; then echo --dbsnp ${inputs["dbsnp"]}; fi) \
		$(if [[ ${inputs["ped"]} != NULL ]]; then echo -ped {inputs["ped"]}; fi) \
		-O ${vcf_dir}/combined.g.vcf \
		|| { echo 'failed to combine gvcf files'; status=1; }

	return $status
}

_genotype_combined_gvcf() {
	local \
		status=0

	$gatk GenotypeGVCFs \
		-R ${inputs["ref"]} \
		-V ${vcf_dir}/combined.g.vcf \
		$(if [[ ${inputs["dbsnp"]} != NULL ]]; then echo --dbsnp ${inputs["dbsnp"]}; fi) \
		$(if [[ ${inputs["ped"]} != NULL ]]; then echo -ped {inputs["ped"]}; fi) \
		-O ${vcf_dir}/genotyped.g.vcf \
		|| { echo 'failed to genotype cohort gvcf file'; status=1; }

	return $status
}

#  emite_gvcf, pemit_gvcf
#
#	* Emite GVCFs from analysis-ready BAM files
#	* I think this subroutine needs renaming...
#	inputs:
#		+ analysis-ready bam files
#
#	tools required:
#		+ gatk HaplotypeCaller
#
#	outputs:
#		+ gvcf files (one per sample)
#
function emit_gvcfs() {
       check_ref; check_gatk_dict; #check_bamlist

       id="$v"
       bam="$(awk '{print $2}' $v)"
       for i in $bam; do
         if [ ! -f "${i/.bam/.bai}" -o ! -f "${i}.bai" ]; then
            echo "Creating SAMTOOLS index: $i"
            samtools index -b -@ $t $i;
         fi
       done   
       #--- Per-sample variant calling emitting GVCFs
       while read -r line; do
            gatk HaplotypeCaller \
            	-R $ref \
            	$(if [[ $ped != NULL ]]; then echo -ped $ped; fi) \
            	$(if [[ $dbsnp != NULL ]]; then echo --dbsnp $dbsnp; fi) \
            	--lenient true \
            	-ERC GVCF \
            	${line}
       done < ${id}
}

function pemit_gvcfs() {
       check_ref; check_gatk_dict; #check_bamlist
       mkdir -p vcall
       id="$v"
       n=$((50/$t))
       bam="$(awk '{print $2}' $v)"
       echo "SAMTOOLS: Index"
       for i in $bam; do
         if [ ! -f "${i/.bam/.bai}" -o ! -f "${i}.bai" ]; then
            echo $i 
         fi
       done > sam.index.list
       cat sam.index.list | parallel --col-sep ' ' echo "index -b -@ $t {}" | xargs -I input -P$n sh -c "samtools input"
       rm sam.index.list
       #--- Per-sample variant calling emitting GVCFs
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo HaplotypeCaller {1} {2}.bam -R $ref $(if [[ $ped != NULL ]]; then echo -ped $ped; fi) $(if [[ $dbsnp != NULL ]]; then echo --dbsnp $dbsnp; fi) -ERC GVCF {3} {4} | xargs -I input -P$n sh -c "gatk input"
}

function vqsr() {
   vcf=$1; op=$2
 
   #SNPs
   gatk VariantRecalibrator \
      -R /mnt/lustre/groups/CBBI1243/KEVIN/db/ucsc.hg19.fasta \
      -V ${vcf} \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/hapmap_3.3.hg19.sites.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/1000G_phase1.indels.hg19.sites.vcf.gz \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/1000G_omni2.5.hg19.sites.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/dbsnp_138.hg19.vcf.gz \
      -an QD \
      -an MQ \
      -an MQRankSum \
      -an ReadPosRankSum \
      -an FS \
      -an SOR \
      -an InbreedingCoeff \
      -mode BOTH \
      -O ${op}.recal \
      --tranches-file ${op}.tranches \
      --rscript-file ${op}.plots.R
   
   #Apply
   gatk ApplyVQSR \
       -V ${vcf} \
       --recal-file ${op}.recal \
       -O ${op}.vqsr-filtered.vcf.gz 

}

#  bcfcall
#
#	* Variant Calling With bcftools
#	inputs:
#		+ list of analysis-ready (marked) bam files
#	tools:
#		+ bcftools (mpileup, index, call, view)
#	outputs:
#		+ a gvcf file containing variants jointly called across the cohort
#
function bcfcall() {
	local \
		tmp_prefix=${tmp_dir}/$(random_id)_ \
		input_bam_stage='marked'

	prep_bamlist $tmp_prefix $input_bam_stage || return 1

	cat ${tmp_prefix}bam.list | awk -v d=${bam_dir}/ '{print d$2}' > ${tmp_prefix}bam.inputs

	$bcftools mpileup \
		--min-MQ 1 \
		--thread ${inputs["threads"]} \
		-f ${inputs["ref"]} \
		-Oz \
		-o ${vcf_dir}/piledup.vcf.gz \
		-b ${tmp_prefix}bam.inputs \
		|| { echo 'bcftools mpileup failed'; return 1; }

	$bcftools index -f -t ${vcf_dir}/piledup.vcf.gz \
		|| { echo 'failed to index bcftools piledup gvcf file'; return 1; }

	$bcftools call \
		-mv \
		--threads ${inputs["threads"]} \
		-Oz \
		-o ${vcf_dir}/called.vcf.gz \
		${vcf_dir}/piledup.vcf.gz \
		|| { echo 'failed to call variants'; return 1; }

	$bcftools index -f -t ${vcf_dir}/called.vcf.gz \
		|| { echo 'failed to index bcftools called gvcf file'; return 1; }

	# unzip the final gvcf file to inspect
	$bcftools view \
		-o ${vcf_dir}/called.vcf \
		${vcf_dir}/called.vcf.gz \

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}


#
#  2. check functions
#
#	summary:
#	+ check_ref()
#	+ check_bwa_idx()
#	+ check_gatk_dict()
#	+ check_samtools_fai()
#	+ check_sites()
#	+ check_adapter()
#	+ check_sample()
#	+ check_fq()
#	+ check_bamlist()
#	+ check_gvcflist()
#

#--- Check References and their indixes
function check_ref() {
       if [[ "${inputs["ref"]}" == NULL ]]; then
          echo -e " -r,--ref not provided! Exiting..."; 1>&2;
          return 1
       elif [ ! -f "${inputs["ref"]}" -o ! -s "${inputs["ref"]}" ]; then
          echo -e " Problem with reference file. Please check that it exists and is not empty..."; 1>&2;
          return 1
       fi
}
function check_bwa_idx() {
       if [[ ! -f "${inputs["ref"]}.bwt" ]]; then
            echo " can't find the bwa index for ${inputs["ref"]}"
			return 1
       fi
}
function check_gatk_dict() {
	if [[ ! -f "${inputs["ref"]/.fasta/.dict}" ]] || [[ ! -f "${inputs["ref"]/.fa/.dict}" ]] ; then
		echo " can't find gatk reference dictionary"
		return 1
	fi
}
function check_samtools_fai() {
       if [[ ! -f "${inputs["ref"]/.fasta/.fai}" ]] || [[ ! -f "${inputs["ref"]/.fa/.fai}" ]]; then
          echo " can't find samtools fai index"
          return 1
       fi
}

#--- Check additional [optional] references (Known sites) for IndelRealignment, BQSR, and VQSR
function check_sites() {
    nks=$(echo $ks | sed 's/,/ --known /g')
    echo "--known $nks" > rtc.ks.txt
    nks=$(echo $ks | sed 's/,/ -known /g')
    echo "-known $nks" > ir.ks.txt
    nks=$(echo $ks | sed 's/,/ --known-sites /g')
    echo "--known-sites $nks" > bqsr.ks.txt
}

#--- Check trimmomatic adapters
function check_adapter() {
    function warning() {
        echo -e """\e[38;5;3mWARNING\e[0m: The adapter was not found! Make sure it is present in the current directory""" 1>&2;
        echo -e """\e[38;5;6m===>\e[0m Attempting to trim without adapter. Press \e[38;5;6mCTRL+C\e[0m to stop\n""" 1>&2;
        return 1;

    }
    case "$(echo "$adap" | tr [:lower:] [:upper:])" in
        NP) if [[ -e "NexteraPE-PE.fa" ]]; then echo ILLUMINACLIP:NexteraPE-PE.fa:2:30:10; else warning; fi ;;
        T3U) if [[ -e "TruSeq3-PE-2.fa" ]]; then echo ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10;  else warning; fi ;;
        T2P) if [[ -e "TruSeq2-PE.fa" ]]; then echo ILLUMINACLIP:TruSeq2-PE.fa:2:30:10;  else warning; fi ;;
        T3P) if [[ -e "TruSeq3-PE.fa" ]]; then echo ILLUMINACLIP:TruSeq3-PE.fa:2:30:10;  else warning; fi ;;
        T2S) if [[ -e "TruSeq2-SE.fa" ]]; then echo ILLUMINACLIP:TruSeq2-SE.fa:2:30:10;  else warning; fi ;;
        T3S) if [[ -e "TruSeq3-SE.fa" ]]; then echo ILLUMINACLIP:TruSeq3-SE.fa:2:30:10;  else warning; fi ;;
         *) echo -e """\e[38;5;3mWARNING\e[0m: No such adapter '$adap'! Type --help for usage\n\e[38;5;6m===>\e[0m Attempting to trim without adapter. Press \e[38;5;6mCTRL+C\e[0m to stop\n""" 1>&2; return 1; ;;
    esac
}

#--- Check sample file
function check_sample() {
	# check meta variabe was set
	if [[ "${inputs["meta"]}" == NULL ]]; then
		echo -e "\e[38;5;1mERROR\e[0m: -s,--sample_list not provided! Exiting..."; 1>&2;
		return 1
	elif [ -f $${inputs["meta"]} -a -s $${inputs["meta"]} ]; then
		for i in $(awk '{print $1}' $${inputs["meta"]}); do
			[ ! $i == "#"* ] && continue
			if [ ! -f ${rds_dir}/$i ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' was not found in the directory '${rds_dir}'.\nPlease specify the path with -p or --path or check that the files in the path are the same in the sample list" 1>&2;
				return 1;
			elif [ -f ${rds_dir}/${i} -a ! -s ${rds_dir}/${i} ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' may be empty. Please check and correct '${rds_dir}'." 1>&2;
				return 1;
		   fi
		done
	elif [ -f ${inputs["meta"]} -a ! -s ${inputs["meta"]} ]; then
		echo -e "\e[38;5;1mERROR\e[0m: '$meta' seems to be empty! Please check and correct." 1>&2;
	fi

	# check read files in $meta file exist
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			[[ -s ${rds_dir}/"${line_array[0]}" ]] || { echo "Error: ${line_array[0]} does not exist or is empty"; return 1; }
			[[ -s ${rds_dir}/"${line_array[1]}" ]] || { echo "Error: ${line_array[0]} does not exist or is empty"; return 1; }
		fi
	done < ${inputs["meta"]}
}

#--- Prepare input for BQSR
function check_gvcflist() {
if [[ ( "$glist" == NULL ) ]]; then
   if [ -f "gvcf.list" ]; then
      rm gvcf.list;
   fi;
   for i in *.gvcf*; do
      if [[ ( -f ${i} ) && ( -s ${i} ) ]]; then  # if gvcf files exist in the current directory and are not empty
         basename -a $(ls $i) >> gvcf.list;
      elif [[ -d vcall ]]; then # if a directory exists called vcall
         for j in vcall/*.gvcf*; do
             if [[ ( -f ${j} ) && ( -s ${j} ) ]]; then # if gvcf files exist in the vcall directory and are not empty
                basename -a $(ls $j) >> gvcf.list;
             fi;
         done;
      else
         echo -e "\n\e[38;5;1mERROR\e[0m: Please check that there are GVCF files in the path $dname\n" 1>&2;
         #return 1;
      fi;
   done;
   if [ -f "gvcf.list" -a -s "gvcf.list" ]; then
      echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat gvcf.list | wc -l) GVCF file(s) counted in '$(readlink -f $(dirname $(cat gvcf.list | head -1)))/' and will be used! Press CTRL+C to stop\n";
      sleep 1;
   fi;
fi
}


#  3. prep functions
#
#	summary:
#	+ prep_fq()
#	+ prep_trim()
#	+ prep_map()
#	+ prep_bamlist()
#	+ prep_gvcflist()
#

#  prep_fq
#
#	* create temporary files used by fq()
#
function prep_fq() { tmp_prefix=$1

	#--- Make input files from forward/reverse runs or SAM/BAM files
	> ${tmp_prefix}fwd.txt
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			echo "${line_array[0]}" >> ${tmp_prefix}fwd.txt
		fi
	done < ${inputs["meta"]}

	> ${tmp_prefix}rev.txt
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			echo "${line_array[1]}" >> ${tmp_prefix}rev.txt
		fi
	done < ${inputs["meta"]}

	if [[ ! -s "${tmp_prefix}rev.txt" ]]; then
		cp ${tmp_prefix}fwd.txt ${tmp_prefix}forward_reverse.txt
		awk -v d="${rds_dir}/" '{print d$1}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}fastq.input.txt
	else
		paste ${tmp_prefix}fwd.txt ${tmp_prefix}rev.txt | awk '{print $1,$2}' > ${tmp_prefix}forward_reverse.txt
		awk -v d="${rds_dir}/" '{print d$1,d$2}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}fastq.input.txt
	fi
	rm ${tmp_prefix}fwd.txt ${tmp_prefix}rev.txt

}

function prep_trim() {
	tmp_prefix=$1
    prep_fq $tmp_prefix || return 1

    #--- On checking for fastq files above, we checked for SAM/BAM as well. If the function picked SAM/BAM, we definitely wanna spill errors since we can't trim SAM/BAM here
    for i in $(awk '{print $1}' ${tmp_prefix}forward_reverse.txt | head -1); do
        if [[ ( ${i} == *.sam ) || ( ${i} == *.sam.gz ) || ( "${i}" == *.bam ) ]]; then
           echo -e "\n\e[38;5;1mERROR\e[0m: No fastq/SAM/BAM file found in the specified location: '${rds_dir}'\nPlease specify path to Fastq/SAM/BAM files using -p or --path\n" 1>&2;
           return 1;
           rm ${tmp_prefix}forward_reverse.txt ${tmp_prefix}fastq.input.txt
        fi
    done
    awk -v d="${rds_dir}/" '{print d$1,d$2, d$1".paired_fp.fq.gz", d$1".unpaired_fu.fq.gz", d$2".paired_rp.fq.gz", d$2".unpaired_ru.fq.gz"}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}trim.input.txt
	#rm forward_reverse.txt fastq.input.txt;
}

#--- Prepare alignment/mapping input
function prep_map() {
	tmp_prefix=$1
	check_ref || return 1
	check_fq $tmp_prefix || return 1
	prep_trim $tmp_prefix || return 1

	if [ -e "${tmp_prefix}forward_reverse.txt" -a -s "${tmp_prefix}forward_reverse.txt" ]; then
		awk -v i="${rds_dir}/" -v o="${sam_dir}/" '{print i$1,i$2,"-o",o$1".sam"}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}align.input.txt
	else
		echo -e "\n\e[38;5;1mERROR\e[0m: Please check that there are fastq files in the path...\n"
	fi
}

#  prep_bamlist
#
#  Prepare input for BQSR, varcall
#	* check that a list of bam files
#	* has been provided by the user
#	* via $blist, 
#	* or create one from the sample list
#	* via $meta
# 
function prep_bamlist() {
	#  * if the user has specified an input bam list 
	#  * then use this one, otherwise, create one
	#  * using the sample IDs in the $inputs["meta"] file.
	local \
		tmp_prefix=$1 \
		input_bam_stage=$2 \
		samples \
		i \
		f

	if [[ ! ${inputs["blist"]}==NULL ]]; then
		echo "using user's input bam list $blist"
		cat ${inputs["blist"]} > ${tmp_prefix}bam.list 
		return 0
	else
		echo "locating bam files attached to the provided sample ids"
		samples=$(sed '/^#/d' ${inputs["meta"]} | awk '{print $3}')
		for i in ${samples[@]}; do
			for f in ${bam_dir}/$i*.${input_bam_stage}.bam; do
				echo "$i $(basename $f)"
			done
		done > ${tmp_prefix}bam.list
	fi

	[[ -s ${tmp_prefix}bam.list ]] || { echo "\n\e[38;5;1mERROR\e[0m: no bam files were found, please check and try again"; return 1; }

}

function prep_gvcflist() {
	#  * if the user has specified an input gvcf list 
	#  * then use this one, otherwise, create one
	#  * using the sample IDs in the $inputs["meta"] file.
	local \
		tmp_prefix=$1 \
		input_gvcf_stage=$2 \
		samples \
		i \
		f

	if [[ ! ${inputs["glist"]}==NULL ]]; then
		echo "using user's input gvcf list $blist"
		cat ${inputs["glist"]} > ${tmp_prefix}gvcf.list 
		return 0
	else
		echo "locating gvcf files attached to the provided sample ids"
		samples=$(sed '/^#/d' ${inputs["meta"]} | awk '{print $3}')
		for i in ${samples[@]}; do
			for f in ${vcf_dir}/$i*.${input_gvcf_stage}.g.vcf; do
				echo "-V $f"
			done
		done > ${tmp_prefix}gvcf.list
	fi

	[[ -s ${tmp_prefix}gvcf.list ]] || { echo "\n\e[38;5;1mERROR\e[0m: no gvcf files were found, please check and try again"; return 1; }

}