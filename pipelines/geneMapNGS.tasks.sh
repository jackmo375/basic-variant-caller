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

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample

	echo ${tmp_prefix}

    check_fq $tmp_prefix
    id=${tmp_prefix}fastq.input.txt; odr="${fqc_dir}/"
    while read -r line; do
        echo -e "running FastQC"
        $fastqc -t $t $line -o $odr
    done < $id

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}

function pfq() {
	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample || return 1
	check_fq $tmp_prefix || return 1

	id=${tmp_prefix}fastq.input.txt; odr="${fqc_dir}/"
	n=$((50/$t))
	echo -e "running FastQC"
	echo -e "Your jobs will be split across $n parallel threads"
	cat $id | \
		parallel --col-sep ' ' echo -e "-t $t {1} $(if [ -n {2} ]; then echo {2}; fi) -o $odr" | \
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
	check_sample || return 1
	prep_trim $tmp_prefix || return 1

	id=${tmp_prefix}trim.input.txt
	while read -r line; do
		echo -e "running Trimmomatic"
		java -jar $trimmomatic PE \
              -phred33 $line \
              $(if [[ $adap != NULL ]]; then checkadapter; fi) \
              LEADING:$leadx \
              TRAILING:$trailx \
              SLIDINGWINDOW:4:15 \
              MINLEN:${minlen} \
              -threads $t
	done < $id
	
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
}

function ptrim() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample || return 1
	prep_trim $tmp_prefix || return 1

	id=${tmp_prefix}trim.input.txt
	n=$((50/$t ))
	echo -e "running Trimmomatic"
	echo -e "Your jobs will be split across $n parallel threads"
	cat $id | \
		parallel --col-sep ' ' echo "PE -phred33 {} $(if [[ $adap != NULL ]]; then checkadapter; fi) LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:$minlen -threads $t" | \
		xargs -I input -P$n sh -c "java -jar $trimmomatic input"
	
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

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample || return 1
	check_ref || return 1
	check_bwa_idx || return 1
	prep_map $tmp_prefix || return 1

	id=${tmp_prefix}align.input.txt
	echo -e "running BWA-BCFTOOLS Alignment/Mapping: Serial\n"

	# align the reads to the reference:
	while read -r line; do
		$bwa mem -t $t $ref $line
	done < $id
	for sam in $(awk '{print $4}' ${tmp_prefix}align.input.txt); do
		# create file names:
		unsorted_bam=${sam/.sam/.bam}
		unsorted_bam=${bam_dir}/${unsorted_bam##*/}
		mapped_bam=${sam/.sam/.mapped.bam}
		mapped_bam=${bam_dir}/${mapped_bam##*/}

		$samtools view \
		   -h \
		   ${sam} \
		   -O BAM \
		   -o $unsorted_bam
		$samtools sort \
		   -O BAM \
		   --reference $ref \
		   -@ $t \
		   -o $mapped_bam \
		   $unsorted_bam

		rm $unsorted_bam
	done

	# remove all temporary files from this call
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*
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
#		+ [gmap()]  ${bam_dir}/<sample_id>.bam (mapped bam files for each sample)
#		+ [pgmap()] ${bam_dir}/<sample_id>_mkdups.bam
#
function gmap() {

	tmp_prefix=${tmp_dir}/$(random_id)_
	check_sample || return 1
	check_ref || return 1
	check_bwa_idx || return 1
	check_gatk_dict || return 1 
	check_samtools_fai || return 1

	id="${tmp_prefix}metadat.txt"
	awk '{print $1,$2,$3,$4}' ${meta} > $id
	echo -e "GATK-BWA Alignment. Your jobs will run in serial\n"

	while read -r line; do
		[[ "$line" == "#"* ]] && continue

		$gatk FastqToSam \
			-F1 ${rds_dir}/$(echo $line | awk '{print $1}') \
			-F2 ${rds_dir}/$(echo $line | awk '{print $2}') \
			-SM $(echo $line | awk '{print $3}') \
			-PL $(echo $line | awk '{print $4}') \
			-RG $(echo $line | awk '{print $3}') \
			-O ${bam_dir}/$(echo $line | awk '{print $3}').unmapped.bam \
			|| { echo 'FastqToSam step failed'; return 1; }
		$bwa mem \
			-t $t \
			$ref $(echo $line | awk -v d="${rds_dir}/" '{print d$1,d$2}') \
			-o ${sam_dir}/$(echo $line | awk '{print d$3}').sam \
			|| { echo 'bwa mem step failed'; return 1; }
		$samtools view \
			-h ${sam_dir}/$(echo $line | awk '{print $3}').sam \
			-O BAM \
			-o ${bam_dir}/$(echo $line | awk '{print $3}').bam \
			|| { echo 'sam to bam step failed'; return 1; }
		$samtools sort \
			-O BAM \
			--reference $ref \
			-@ $t \
			-o ${bam_dir}/$(echo $line | awk '{print $3}').mapped.bam \
			${bam_dir}/$(echo $line | awk '{print $3}').bam \
			|| { echo 'sorting the bam files step failed'; return 1; }
		$gatk MergeBamAlignment \
			-O ${bam_dir}/$(echo $line | awk '{print $3}').bam \
			-R ${ref} \
			-UNMAPPED ${bam_dir}/$(echo $line | awk '{print $3}').unmapped.bam \
			-ALIGNED ${bam_dir}/$(echo $line | awk '{print $3}').mapped.bam \
			|| { echo 'merge bam alignment step failed'; return 1; }

		for i in ${bam_dir}/$(echo $line | awk '{print $3}').*mapped.bam; do
			if [ -e $i ]; then
				rm $i
			fi
		done

		# do we not need to mark duplicates here too?

	done < $id

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

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
#		+ bam list
function bqsr() {

	check_ref || return 1
	check_gatk_dict || return 1
	check_bamlist || return 1

	id="$blist"
	n=$((50/$t))

	while read -r line; do
		$gatk BaseRecalibrator \
			-I ${dname}/${line/.bam/} \
			-R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "bqsr.ks.txt" -a -s "bqsr.ks.txt" ]; then rm rtc.ks.txt ir.ks.txt; cat bqsr.ks.txt; fi; fi) \
			-O bqsr/${line/.bam/_recal_data.table}
		$gatk ApplyBQSR \
			-R $ref 
			-I ${dname}/${line/.bam/} 
			--bqsr-recal-file bqsr/${line/.bam/_recal_data.table} \
			-O aligned/${line/.bam/.mapped.bam}
	done < ${id}
	if [ -e "${id}" ]; then rm ${id}; fi
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

#--- Emite GVCFs from analysis-ready BAM files
function emit_gvcfs() {
       check_ref; check_gatk_dict; #check_bamlist
       mkdir -p vcall
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
            gatk HaplotypeCaller -R $ref $(if [[ $ped != NULL ]]; then echo -ped $ped; fi) $(if [[ $dbsnp != NULL ]]; then echo --dbsnp $dbsnp; fi) --lenient true -ERC GVCF ${line}
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

function combinegvcfs() {
       gatk CombineGVCFs \
          -R $ref \
          $(if [[ $dbsnp != NULL ]]; then echo --dbsnp $dbsnp; fi) \
          $(if [[ $ped != NULL ]]; then echo -ped $ped; fi) \
          --arguments_file ${v} \
          -O vcall/${out}.gvcf.gz 
}

function genogvcfs() {
       gatk GenotypeGVCFs \
          -R $ref \
          $(if [[ $dbsnp != NULL ]]; then echo --dbsnp $dbsnp; fi) \
          $(if [[ $ped != NULL ]]; then echo -ped $ped; fi) \
          --arguments_file ${v} \
          -O vcall/${out}.vcf.gz 
}

#--- Variant Calling with GATK (Single Cohort Joint) in Serial
function varcall() {
         if [[ ( $glist == NULL ) && ( $blist == NULL ) ]]; then
            echo -e "\e[38;5;3mWARNING\e[0m: Neither -g,--gvcf_list nor -b,--bam_list provided... "
            if [ -f "gvcf.list" ]; then rm gvcf.list; elif [ -f "bam.list" ]; then rm bam.list; fi
            sleep 1;
            if [ -d "vcall" -a -n "$(ls -A vcall)" ]; then 
               echo -e "Checking GVCF files in $(readlink -f vcall)/..."; 
               sleep 1;
               for j in vcall/*.gvcf*; do
                   if [ -f ${j} -a -s ${j} ]; then # if gvcf files exist in the vcall directory and are not empty
                      echo "-V $j";
                   fi
               done | sed '/.tbi/d' > gvcf.list
               dname="$(readlink -f vcall)/"
               if [ -f "gvcf.list" -a -s "gvcf.list" ]; then
                  echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat gvcf.list | wc -l) valid GVCF file(s) counted in '$(readlink -f $dname)/' and will be used! Press CTRL+C to stop\n";
                  sleep 1;
                  v="gvcf.list"
                 #check index of bam files
                 gv="$(awk '{print $2}' $v)"
                 for i in ${gv}; do
                   if [[ ! -f "${i}.tbi" ]]; then
                      echo "Generating VCF index: $i" 
                      tabix -p vcf -f $i;
                   fi
                 done
                  check_ref
		  if [[ $(cat $v | wc -l) != 1 ]]; then
                     combinegvcfs;
                     rm $v
		     echo "-V vcall/${out}.gvcf.gz" > genogvcf.in.txt
		     v="genogvcf.in.txt"
		  fi
                  genogvcfs;
                  rm $v
               fi
            elif [[ -d "aligned" ]]; then
               echo -e "Checking BAM files in $(readlink -f aligned)/...";
               sleep 1;
               for j in aligned/*.bam; do
                   if [ -f ${j} -a -s ${j} ]; then # if gvcf files exist in the vcall directory and are not empty
                      echo "-I $j -O vcall/$(basename ${j/.bam/.gvcf.gz})";
                   fi
               done > bam.list
               dname="$(readlink -f aligned)/"
               if [ -f "bam.list" -a -s "bam.list" ]; then
                  echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat bam.list | wc -l) valid BAM file(s) counted in '$(readlink -f $dname)/' and will be used! Press CTRL+C to stop\n";
                  sleep 1;
                  v="bam.list"
                  check_ref
                  if [[ $cmd == "gatkcall" ]]; then
                      emit_gvcfs; 
                      awk '{print "-V",$4}' $v > gvcf.list; v="gvcf.list"; #Make GVCF list from BAM list after emmitting GVCFs to pass to CombineGVCFs
                      combinegvcfs; genogvcfs
                      rm $v
                  elif [[ $cmd == "pgatkcall" ]]; then
                      pemit_gvcfs; 
                      awk '{print "-V",$4}' $v > gvcf.list; v="gvcf.list"; #Make GVCF list from BAM list after emmitting GVCFs to pass to CombineGVCFs
                      combinegvcfs; genogvcfs
                      rm $v
                  fi
               fi
            fi
         elif [ -f $glist -a -s $glist ]; then
              mkdir -p vcall
              #Get directory name of the gvcf files
              dn="$(readlink -f $(dirname $(cat $glist | head -1)))/"
              if [[ "$dname" == NULL ]]; then 
                 dname="$dn"; 
              fi

             #Get basename of all gvcf files
             for i in $(cat $glist); do
                 if [ -f ${dname}$(basename ${i}) -a -s ${dname}$(basename ${i}) ]; then
                    echo $(basename $i);
                 fi
              done > gb.list

              for i in $(cat gb.list); do
                 if [ -f ${dname}${i} -a -s ${dname}${i} ]; then
                    echo "-V ${dname}${i}";
                 fi
              done > g.list
              mv g.list ${glist/.*/.gvcfs.in.txt}
              v="${glist/.*/.gvcfs.in.txt}"
              if [ -e gb.list ]; then rm gb.list; fi
              if [ -f "$v" -a -s "$v" ]; then
                 echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat $glist | wc -l) valid GVCF file(s) counted in '$(readlink -f $dname)/' and will be used! Press CTRL+C to stop\n";
                 sleep 1;
                 check_ref
                 for i in $(cut -f2 -d' ' ${v}); do
                   if [[ ! -f "${i}.tbi" ]]; then
                      echo "Generating VCF index: $i" 
                      tabix -p vcf -f $i;
                   fi
                 done
                 if [[ $(cat $v | wc -l) != 1 ]]; then
                     combinegvcfs;
                     rm $v
                     echo "-V vcall/${out}.gvcf.gz" > genogvcf.in.txt
                     v="genogvcf.in.txt"
                  fi
                  genogvcfs;
                  rm $v
              elif [ -f "$v" -a ! -s "$v" ]; then
                 rm $v
                 echo -e "\e[38;5;1mERROR\e[0m: Problem with gvcf list. Did you forget to specify the [CORRECT] path to gvcf file(s) with -p,--path ?" 1>&2;
                 return 1;
              fi
              if [ -e $v ]; then rm $v; fi
         elif [ -f $blist -a -s $blist ]; then

              #Get directory name of the bam files
              dn="$(readlink -f $(dirname $(cat $blist | head -1)))/"
              if [[ "$dname" == NULL ]]; then 
                 dname="$dn"; 
              fi

             #Get basename of all bam files
              for i in $(cat $blist); do
                 if [ -f  ${dname}$(basename ${i}) -a -s  ${dname}$(basename ${i}) ]; then
                    echo $(basename $i);
                 fi
              done > bb.list

              for i in $(cat bb.list); do
                 if [ -f ${dname}${i} -a -s ${dname}${i} ]; then
                    echo "-I ${dname}${i} -O vcall/$(basename ${i/.bam/.gvcf.gz})"
                 fi
              done > b.list
              mv b.list ${blist/.*/.hapcaller.in.txt}
              v="${blist/.*/.hapcaller.in.txt}"
              if [ -f bb.list ]; then rm bb.list; fi
              if [ -f "$v" -a -s "$v" ]; then
                 echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat $blist | wc -l) valid BAM file(s) counted in '$(readlink -f $dname)/' and will be used! Press CTRL+C to stop\n";
                 sleep 1;
                  check_ref
                  if [[ $cmd == "gatkcall" ]]; then
                      emit_gvcfs; 
                      awk '{print "-V",$4}' $v > gvcf.list; v="gvcf.list"; #Make GVCF list from BAM list after emmitting GVCFs to pass to CombineGVCFs
                      combinegvcfs; 
		      genogvcfs
                  elif [[ $cmd == "pgatkcall" ]]; then
                      pemit_gvcfs; 
                      awk '{print "-V",$4}' $v > gvcf.list; v="gvcf.list"; #Make GVCF list from BAM list after emmitting GVCFs to pass to CombineGVCFs
                      combinegvcfs; 
		      genogvcfs
                  fi
                 if [ -f $v ]; then rm $v; fi
              elif [ -f "$v" -a ! -s "$v" ]; then
                 rm $v
                 echo -e "\e[38;5;1mERROR\e[0m: Problem with bam list. Did you forget to specify the [CORRECT] path to bam file(s) with -p,--path ?" 1>&2;
                 return 1;
              fi
         elif [[ ( ( -f $glist ) && ( ! -s $glist ) ) || ( ( -f $blist ) && ( ! -s $blist ) ) ]]; then
              echo -e "\e[38;5;1mERROR\e[0m: Problem with [gvcf/bam] list. Please check and correct." 1>&2;
              return 1;
         else       
             gcallhelp 1>&2;
             return 1;
         fi
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
#--- Variant Calling With bcftools
function bcfcall() {
       check_ref; check_bamlist
       mkdir -p vcall
       id="$blist"

       dn="$(readlink -f $(dirname $(cat $blist | head -1)))/"
       if [[ "$dname" == NULL ]]; then
          dname="$dn";
       fi
       #Get basename of all bam files
       for i in $(cat $blist); do
          if [ -f  ${dname}$(basename ${i}) -a -s  ${dname}$(basename ${i}) ]; then
             echo $(basename $i);
          fi
       done > bb.list
       #attach directory name to bam files
       for i in $(cat bb.list); do
          if [ -f "${dname}${i}" -a -s "${dname}${i}" ]; then
             echo "${dname}$(basename ${i})"
          fi
       done > b.list
       v="b.list"

      if [ -f "$v" -a -s "$v" ]; then
         echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat $blist | wc -l) valid BAM file(s) counted in '$(readlink -f $dname)/' and will be used! Press CTRL+C to stop\n";
         sleep 1;
          check_ref;
          bcftools mpileup --min-MQ 1 --thread $t -f $ref -Oz -o out.vcf.gz -b $v
          bcftools index -f -t out.vcf.gz
          bcftools call -mv --threads $t -Oz -o vcall/${out}.vcf.gz out.vcf.gz
          bcftools index -f -t vcall/${out}.vcf.gz
          for i in $v out.vcf.gz; do if [[ -e "${i}" ]]; then rm $i; fi; done
      elif [ -f "$v" -a ! -s "$v" ]; then
         rm $v
         echo -e "\e[38;5;1mERROR\e[0m: Problem with bam list. Did you forget to specify the [CORRECT] path to bam file(s) with -p,--path ?" 1>&2;
         return 1;
      fi

#      echo -e "Variant Calling - BCFTOOLS"
#      bcftools mpileup --min-MQ 1 --thread $t -f $ref -Oz -o out.vcf.gz -b $id
#      bcftools index -f -t out.vcf.gz
#      bcftools call -mv --threads $t -Oz -o vcall/${out}.vcf.gz out.vcf.gz
#      bcftools index -f -t vcall/${out}.vcf.gz
#      for i in $id out.vcf.gz; do if [[ -e "${i}" ]]; then rm $i; fi; done
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
       if [[ "$ref" == NULL ]]; then
          echo -e "\e[38;5;1mERROR\e[0m: -r,--ref not provided! Exiting..."; 1>&2;
          return 1
       elif [ ! -f "$ref" -o ! -s "$ref" ]; then
          echo -e "\e[38;5;1mERROR\e[0m: Problem with reference file. Please check that it exists and is not empty..."; 1>&2;
          return 1
       fi
}
function check_bwa_idx() {
       check_ref
       if [[ ! -f "${ref}.bwt" ]]; then
            echo "${error}: can't find the bwa index for $ref"
			return 1
       fi
}
function check_gatk_dict() {
	check_ref
	if [[ ! -f "${ref/.fasta/.dict}" ]]; then
		echo "${error}: can't find gatk reference dictionary"
		return 1
	fi
}
function check_samtools_fai() {
       check_ref
       if [[ ! -f "${ref/.fasta/.fai}" ]]; then
          echo "${error}: can't find samtools fai index"
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
	if [[ "$meta" == NULL ]]; then
		echo -e "\e[38;5;1mERROR\e[0m: -s,--sample_list not provided! Exiting..."; 1>&2;
		return 1
	elif [ -f $meta -a -s $meta ]; then
		for i in $(awk '{print $1}' $meta); do
			[ ! $i == "#"* ] && continue
			if [ ! -f ${rds_dir}$i ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' was not found in the directory '${rds_dir}'.\nPlease specify the path with -p or --path or check that the files in the path are the same in the sample list" 1>&2;
				return 1;
			elif [ -f ${rds_dir}${i} -a ! -s ${rds_dir}${i} ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' may be empty. Please check and correct '${rds_dir}'." 1>&2;
				return 1;
		   fi
		done
	elif [ -f $meta -a ! -s $meta ]; then
		echo -e "\e[38;5;1mERROR\e[0m: '$meta' seems to be empty! Please check and correct." 1>&2;
	fi
}

#--- Check Fastq/SAM/BAM and make input files
function check_fq() {
	tmp_prefix=$1

	# get all sample fastq file names from meta file

    #--- Make input files from forward/reverse runs or SAM/BAM files
	> ${tmp_prefix}fwd.txt
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
    		echo "${line_array[0]}" >> ${tmp_prefix}fwd.txt
		fi
	done < $meta
    if [[ ! -s "${tmp_prefix}fwd.txt" ]]; then
       echo -e "\n\e[38;5;1mERROR\e[0m: No fastq/SAM/BAM file found in the specified location: '$dname'\nPlease specify path to Fastq/SAM/BAM files using -p or --path\n"
       rm ${tmp_prefix}fwd.txt 1>&2;
       return 1;
    fi

	> ${tmp_prefix}rev.txt
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
    		echo "${line_array[1]}" >> ${tmp_prefix}rev.txt
		fi
	done < $meta

	if [[ ! -s "${tmp_prefix}rev.txt" ]]; then
		cp ${tmp_prefix}fwd.txt ${tmp_prefix}forward_reverse.txt
		awk -v d="${rds_dir}/" '{print d$1}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}fastq.input.txt
	else
		paste ${tmp_prefix}fwd.txt ${tmp_prefix}rev.txt | awk '{print $1,$2}' > ${tmp_prefix}forward_reverse.txt
		awk -v d="${rds_dir}/" '{print d$1,d$2}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}fastq.input.txt
	fi
    rm ${tmp_prefix}fwd.txt ${tmp_prefix}rev.txt

}

#--- Prepare input for BQSR
function check_bamlist() {
	tmp_prefix=$1

	[[ ! "$blist" == NULL ]] && return 0

	for i in *.bam; do
	  if [[ ( -f ${i} ) && ( -s ${i} ) ]]; then  # if bam files exist in the current directory and are not empty
		 basename -a $(ls $i) >> bam.list;
	  elif [[ -d aligned ]]; then # if a directory exists called aligned
		 for j in aligned/*.bam; do
		     if [[ ( -f ${j} ) && ( -s ${j} ) ]]; then # if bam files exist in the aligned directory and are not empty
		        basename -a $(ls $j) >> bam.list;
		     fi;
		 done;
	  else
		 echo -e "\n\e[38;5;1mERROR\e[0m: Please check that there are bam files in the path $dname\n" 1>&2;
	 return 1;
	  fi;
	done;
	if [ -f "bam.list" -a -s "bam.list" ]; then
	  echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat bam.list | wc -l) BAM file(s) counted in '$(readlink -f $(dirname $(cat bam.list | head -1)))/' and will be used! Press CTRL+C to stop\n";
	  sleep 1;
	fi;
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
#	+ prep_trim()
#	+ prep_map()
#


function prep_trim() {
	tmp_prefix=$1
    check_fq $tmp_prefix || return 1

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

