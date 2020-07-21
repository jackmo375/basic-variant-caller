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
###############################################

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

#--- FastQC
function fq() {
    checkfq
    id=fastq.input.txt; odr="${fqc_dir}/"
    while read -r line; do
        echo -e "FastQC"
        $fastqc -t $t $line -o $odr
    done < $id
    rm fastq.input.txt
}

function pfq() {
       checkfq
       id=fastq.input.txt; odr="${fqc_dir}/"
       n=$((50/$t))
       echo -e "FastQC"
       echo -e "Your jobs will be run in $n parallel runs"
       cat $id | \
           parallel --col-sep ' ' echo -e "-t $t {1} $(if [ -n {2} ]; then echo {2}; fi) -o $odr" | \
           xargs -I input -P$n sh -c "$fastqc input"
       rm fastq.input.txt
}

#--- trimmomatic
function trim() {
       preptrim
       mkdir -p paired unpaired
       id=trim.input.txt
       while read -r line; do
           echo -e "Trommomatic"
           trimmomatic PE \
              -phred33 $line \
              $(if [[ $adap != NULL ]]; then checkadapter; fi) \
              LEADING:$leadx \
              TRAILING:$trailx \
              SLIDINGWINDOW:4:15 \
              MINLEN:36 \
              -threads $t
       done < $id
       rm trim.input.txt
}

function ptrim() {
       preptrim
       mkdir -p paired unpaired
       id=trim.input.txt
       n=$((50/$t ))
       echo -e "trimmomatic"
       echo -e "Your jobs will be run in $n parallel runs"
          cat $id | \
              parallel --col-sep ' ' echo "PE -phred33 {} $(if [[ $adap != NULL ]]; then checkadapter; fi) LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:36 -threads $t" | \
              xargs -I input -P$n sh -c "trimmomatic input"
       rm trim.input.txt
}

#--- Mapping/Alignment (BWA)
function bmap() {
       prepmap
       mkdir -p aligned
       id=align.input.txt
       echo -e "BWA-BCFTOOLS Alignment/Mapping: Serial\n"
       while read -r line; do
           $bwa mem -t $t ${ref_dir}/${species}.bwa $line
       done < $id
       for sam in $(awk '{print $4}' align.input.txt); do
           $samtools view \
               -h \
               ${sam} \
               -O BAM \
               -o ${sam/.sam/.bam}
           $samtools sort \
               -O BAM \
               --reference $ref \
               -@ $t \
               -o ${sam/.sam/.mapped.bam} \
               ${sam/.sam/.bam}
           echo ${sam/.sam/.mapped.bam}
           rm ${sam/.sam/.bam}
       done > bam.list
       rm align.input.txt
}

function pbmap() {
       prepmap
       mkdir -p aligned
       id=align.input.txt
       awk '{print $4,"-O BAM -o",$4}' align.input.txt | \
           sed 's/.sam/.bam/2' > sam2bam.input.txt
       awk '{print $5,$5}' sam2bam.input.txt | \
           sed 's/.bam/.mapped.bam/1' > sortbam.input.txt
       n=$((50/$t))
       echo -e "BWA-BCFTOOLS Alignment/Mapping\n"
       cat align.input.txt | \
           parallel --col-sep ' ' echo "mem -t $t ${ref_dir}/${species}.bwa {}" | \
           xargs -I input -P$n sh -c "$bwa input"
       cat sam2bam.input.txt | \
           parallel --col-sep ' ' echo "view -h {}" | \
           xargs -I input -P$n sh -c "$samtools input"
       cat sortbam.input.txt | \
           parallel --col-sep ' ' echo "sort -O BAM --reference $ref -@ $t -o {}" | \
           xargs -I input -P$n sh -c "$samtools input"
       for sam in $(awk '{print $4}' align.input.txt); do
           rm ${sam};
       done
       for bam in $(awk '{print $2}' sortbam.input.txt); do
           rm ${bam};
       done
       for i in out.vcf.gz aligned/*.sam; do
           if [[ -e "${i}" ]]; then
              rm $i;
           fi;
       done
       rm align.input.txt sam2bam.input.txt sortbam.input.txt
}

#--- GATKv4 BWA Mapping/Alignment
function gmap() {
       check_sample; check_ref; check_bwa_idx; check_gatk_dict; check_samtools_fai;
       mkdir -p aligned
       awk '{print $1,$2,$3,$4}' ${meta} > metadat.txt
       id="metadat.txt"
       echo -e "GATK-BWA Alignment. Your jobs will run in serial\n"
       while read -r line; do
           $gatk FastqToSam \
                -F1 ${dname}$(echo $line | awk '{print $1}') \
                -F2 ${dname}$(echo $line | awk '{print $2}') \
                -SM $(echo $line | awk '{print $3}') \
                -PL $(echo $line | awk '{print $4}') \
                -RG $(echo $line | awk '{print $3}') \
                -O aligned/"$(echo $line | awk '{print $3}').unmapped.bam"
           $bwa mem \
                -t $t \
                $ref $(echo $line | awk '{print $1,$2}') \
                -o "$(echo $line | awk '{print "aligned/"$3}').sam"
           $samtools view \
                -h "$(echo $line | awk '{print "aligned/"$3}').sam" \
                -O BAM \
                -o "$(echo $line | awk '{print "aligned/"$3}').bam"
           $samtools sort \
                -O BAM \
                --reference $ref \
                -@ $t \
                -o "$(echo $line | awk '{print "aligned/"$3}').mapped.bam" \
                "$(echo $line | awk '{print "aligned/"$3}').bam"
           $gatk MergeBamAlignment \
                -O "$(echo $line | awk '{print "aligned/"$3}').bam" \
                -R ${ref} \
                -UNMAPPED "$(echo $line | awk '{print "aligned/"$3}').unmapped.bam" \
                -ALIGNED "$(echo $line | awk '{print "aligned/"$3}').mapped.bam"
       done < $id
       for i in $(echo $line | awk '{print "aligned/"$3}').*mapped.bam; do
           if [ -e $i ]; then
              rm $i
           fi
       done
       rm metadat.txt
}

function pgmap() {
       check_sample; check_ref; check_bwa_idx; check_gatk_dict; check_samtools_fai;
       mkdir -p aligned
       awk '{print $1,$2,$3,$4}' ${meta} > metadat.txt
       id="metadat.txt"
       n=$((50/$t))
       echo -e "GATK-BWA Alignment and Mark Duplicates. Your jobs will be run in $n parallel runs\n"
       #--- Make unmapped BAM files from raw FASTQ files
       cat ${id} | parallel --col-sep ' ' echo "FastqToSam -F1 ${dname}{1} -F2 ${dname}{2} -SM {3} -PL {4} -RG {3} -O aligned/{3}.unmapped.bam" | xargs -I input -P$n sh -c "$gatk input"
       cat ${id} | parallel --col-sep ' ' echo "mem -t $t $ref ${dname}{1} ${dname}{2} -o aligned/{3}.mapped.sam" | xargs -I input -P$n sh -c "$bwa input"
       cat ${id} | parallel --col-sep ' ' echo "view -O BAM -h aligned/{3}.mapped.sam -o aligned/{3}.unsorted.mapped.bam" | xargs -I input -P$n sh -c "$samtools input"
       cat ${id} | parallel --col-sep ' ' echo "sort -O BAM --reference $ref -@ $t -o aligned/{3}.mapped.bam aligned/{3}.unsorted.mapped.bam" | xargs -I input -P$n sh -c "$samtools input"
       rm aligned/*.sam aligned/*.unsorted.mapped.bam
       cat ${id} | parallel --col-sep ' ' echo "MergeBamAlignment -O aligned/{3}.bam -R ${ref} -UNMAPPED aligned/{3}.unmapped.bam -ALIGNED aligned/{3}.mapped.bam" | xargs -I input -P$n sh -c "$gatk input"
       rm aligned/*mapped.bam
       #--- Mark duplicates (and remove)
       #mkdir -p mkd
       cat ${id} | parallel --col-sep ' ' echo MarkDuplicates -I aligned/{3}.bam -O aligned/{3}_mkdups.bam -M aligned/{3}_marked_dup_metrics.txt --REMOVE_DUPLICATES false | xargs -I input -P$n sh -c "$gatk input"
       cat ${id} | parallel --col-sep ' ' echo index -b aligned/{3}_mkdups.bam | xargs -I input -P$n sh -c "$samtools input"
       cat ${id} | parallel --col-sep ' ' rm aligned/{3}.bam
}

#--- Indel Realignment (THis will not be run if GATK is used since HaplotypeCaller essentially does local rearrangements)

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

#--- Base Quality Score Recallibration (BQSR)
function bqsr() {
       if [[ "$dname" == NULL ]]; then
          echo -e "\n\e[38;5;1mERROR\e[0m: -p,--path not provided! Please specify path to BAM files"; 1>&2;
          exit 1;
       fi
       check_ref; check_gatk_dict; check_bamlist
       id="$blist"
       n=$((50/$t))
       mkdir -p bqsr; mkdir -p aligned
       while read -r line; do
             gatk BaseRecalibrator -I ${dname}/${line/.bam/} -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "bqsr.ks.txt" -a -s "bqsr.ks.txt" ]; then rm rtc.ks.txt ir.ks.txt; cat bqsr.ks.txt; fi; fi) -O bqsr/${line/.bam/_recal_data.table}
             gatk ApplyBQSR -R $ref -I ${dname}/${line/.bam/} --bqsr-recal-file bqsr/${line/.bam/_recal_data.table} -O aligned/${line/.bam/.mapped.bam}
       done < ${id}
       if [ -e "${id}" ]; then rm ${id}; fi
}

function pbqsr() {
       if [[ "$dname" == NULL ]]; then
          echo -e "\n\e[38;5;1mERROR\e[0m: -p,--path not provided! Please specify path to BAM files"; 1>&2;
          exit 1;
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
                 exit 1;
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
                 exit 1;
              fi
         elif [[ ( ( -f $glist ) && ( ! -s $glist ) ) || ( ( -f $blist ) && ( ! -s $blist ) ) ]]; then
              echo -e "\e[38;5;1mERROR\e[0m: Problem with [gvcf/bam] list. Please check and correct." 1>&2;
              exit 1;
         else       
             gcallhelp 1>&2;
             exit 1;
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
         exit 1;
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
#	+ checkfq()
#	+ check_bamlist()
#	+ check_gvcflist()
#

#--- Check References and their indixes
function check_ref() {
       if [[ "$ref" == NULL ]]; then
          echo -e "\e[38;5;1mERROR\e[0m: -r,--ref not provided! Exiting..."; 1>&2;
          exit 1
       elif [ ! -f "$ref" -o ! -s "$ref" ]; then
          echo -e "\e[38;5;1mERROR\e[0m: Problem with reference file. Please check that it exists and is not empty..."; 1>&2;
          exit 1
       fi
}
function check_bwa_idx() {
       check_ref
       if [[ ! -f "${ref}.bwt" ]]; then
            $bwa index $ref
       fi
}
function check_gatk_dict() {
       check_ref
       #--- Create a reference dictionary if it does not exist
       if [[ ! -f "${ref/.fasta/.dict}" ]]; then
          $gatk CreateSequenceDictionary -R $ref
       fi
}
function check_samtools_fai() {
       check_ref
       if [[ ! -f "${ref/.fasta/.fai}" ]]; then
          $samtools faidx $ref
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
function checkadapter() {
    function warning() {
        echo -e """\e[38;5;3mWARNING\e[0m: The adapter was not found! Make sure it is present in the current directory""" 1>&2;
        echo -e """\e[38;5;6m===>\e[0m Attempting to trim without adapter. Press \e[38;5;6mCTRL+C\e[0m to stop\n""" 1>&2;
        exit 1;

    }
    case "$(echo "$adap" | tr [:lower:] [:upper:])" in
        NP) if [[ -e "NexteraPE-PE.fa" ]]; then echo ILLUMINACLIP:NexteraPE-PE.fa:2:30:10; else warning; fi ;;
        T3U) if [[ -e "TruSeq3-PE-2.fa" ]]; then echo ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10;  else warning; fi ;;
        T2P) if [[ -e "TruSeq2-PE.fa" ]]; then echo ILLUMINACLIP:TruSeq2-PE.fa:2:30:10;  else warning; fi ;;
        T3P) if [[ -e "TruSeq3-PE.fa" ]]; then echo ILLUMINACLIP:TruSeq3-PE.fa:2:30:10;  else warning; fi ;;
        T2S) if [[ -e "TruSeq2-SE.fa" ]]; then echo ILLUMINACLIP:TruSeq2-SE.fa:2:30:10;  else warning; fi ;;
        T3S) if [[ -e "TruSeq3-SE.fa" ]]; then echo ILLUMINACLIP:TruSeq3-SE.fa:2:30:10;  else warning; fi ;;
         *) echo -e """\e[38;5;3mWARNING\e[0m: No such adapter '$adap'! Type --help for usage\n\e[38;5;6m===>\e[0m Attempting to trim without adapter. Press \e[38;5;6mCTRL+C\e[0m to stop\n""" 1>&2; exit 1; ;;
    esac
}

#--- Check sample file
function check_sample() {
   if [[ "$meta" == NULL ]]; then
      echo -e "\e[38;5;1mERROR\e[0m: -s,--sample_list not provided! Exiting..."; 1>&2;
      exit 1
   elif [ -f $meta -a -s $meta ]; then
        for i in $(awk '{print $1}' $meta); do
           if [ ! -f ${dname}$i ]; then
	      echo -e "\e[38;5;1mERROR\e[0m: '${i}' was not found in the directory '${dname}'.\nPlease specify the path with -p or --path or check that the files in the path are the same in the sample list" 1>&2;
              exit 1;
           elif [ -f ${dname}${i} -a ! -s ${dname}${i} ]; then
              echo -e "\e[38;5;1mERROR\e[0m: '${i}' may be empty. Please check and correct '${dname}'." 1>&2;
              exit 1;
           fi
        done
   elif [ -f $meta -a ! -s $meta ]; then
        echo -e "\e[38;5;1mERROR\e[0m: '$meta' seems to be empty! Please check and correct." 1>&2;
   fi
}

#--- Check Fastq/SAM/BAM and make input files
function checkfq() {
    #--- Make input files from forward/reverse runs or SAM/BAM files
    for i in ${dname}*_1.fastq* ${dname}*_R1*.fastq* ${dname}*_1.fq* ${dname}*_R1*.fq* ${dname}*.1.fq* ${dname}*.R1.fq* ${dname}*.sam ${dname}*.bam; do
        if [[ -e $i ]]; then
           basename $i;
        fi
    done > fwd.txt
    if [[ ! -s "fwd.txt" ]]; then
       echo -e "\n\e[38;5;1mERROR\e[0m: No fastq/SAM/BAM file found in the specified location: '$dname'\nPlease specify path to Fastq/SAM/BAM files using -p or --path\n"
       rm fwd.txt 1>&2;
       exit 1;
    fi

    for i in ${dname}*_2.fastq* ${dname}*_R2*.fastq* ${dname}*_2.fq* ${dname}*_R2*.fq* ${dname}*.2.fq* ${dname}*.R2.fq*; do
        if [[ -e $i ]]; then
           basename $i;
        fi
    done > rev.txt
     if [[ ! -s "rev.txt" ]]; then
        cp fwd.txt forward_reverse.txt
        awk -v d="${dname}" '{print d$1}' forward_reverse.txt > fastq.input.txt
     else
        paste fwd.txt rev.txt | awk '{print $1,$2}' > forward_reverse.txt
        awk -v d="${dname}" '{print d$1,d$2}' forward_reverse.txt > fastq.input.txt
     fi
    rm fwd.txt rev.txt
}

#--- Prepare input for BQSR
function check_bamlist() {
if [[ "$blist" == NULL ]]; then
   if [ -f "bam.list" ]; then
      rm bam.list;
   fi;
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
	 exit 1;
      fi;
   done;
   if [ -f "bam.list" -a -s "bam.list" ]; then
      echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat bam.list | wc -l) BAM file(s) counted in '$(readlink -f $(dirname $(cat bam.list | head -1)))/' and will be used! Press CTRL+C to stop\n";
      sleep 1;
   fi;
fi
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
         #exit 1;
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
#	+ preptrim()
#	+ prepmap()
#

#--- Check if Fastq files are present and make trim input files
function preptrim() {
    checkfq
    #--- On checking for fastq files above, we checked for SAM/BAM as well. If the function picked SAM/BAM, we definitely wanna spill errors since we can't trim SAM/BAM here
    for i in $(awk '{print $1}' forward_reverse.txt | head -1); do
        if [[ ( ${i} == *.sam ) || ( ${i} == *.sam.gz ) || ( "${i}" == *.bam ) ]]; then
           echo -e "\n\e[38;5;1mERROR\e[0m: No fastq/SAM/BAM file found in the specified location: '$dname'\nPlease specify path to Fastq/SAM/BAM files using -p or --path\n" 1>&2;
           exit 1;
           rm forward_reverse.txt fastq.input.txt
        fi
    done
    mkdir -p paired
    awk -v d="${dname}" '{print d$1,d$2,"paired/"$1"_fp.fq.gz","unpaired/"$1"_fu.fq.gz","paired/"$2"_rp.fq.gz","unpaired/"$2"_ru.fq.gz"}' forward_reverse.txt > trim.input.txt
    #rm forward_reverse.txt fastq.input.txt;
}

#--- Prepare alignment/mapping input
function prepmap() {
   check_ref; checkfq; preptrim
   if [ -e "forward_reverse.txt" -a -s "forward_reverse.txt" ]; then
      awk -v d="$dname" '{print d$1,d$2,"-o","aligned/"$1".sam"}' forward_reverse.txt > align.input.txt
      #rm trim.input.txt
   else
      echo -e "\n\e[38;5;1mERROR\e[0m: Please check that there are fastq files in the path...\n"
   fi
}

