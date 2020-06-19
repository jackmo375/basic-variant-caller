#!/bin/bash

source ../configs/pipeline.cfg
source ../configs/settings.cfg

if [ $# -lt 5 ];then
	echo "Usage: ./alignment_single.sh [aligner] [reference ver.] [prefix] [fastq1] [fastq2]"
fi

sw_id=$1 # argument 1: the name of aligner
ref_id=$2 # argument 2: the version of reference

sample_id=$3 # argument 3: prefix for outputs

fastq1=$4 # argument 4,5: input fastq files (gzipped)
fastq2=$5

prefix="${sample_id}_${sw_id}_${ref_id}"

if [ ! -d  ${tmp_dir}/${prefix} ];then
	mkdir ${tmp_dir}/${prefix}
fi

if [ ! -d ${work_dir}/${prefix} ];then
	mkdir ${work_dir}/${prefix}
	mkdir ${work_dir}/${prefix}/log
fi

LOG_FILE="${work_dir}/${prefix}/log/${prefix}.log"

echo "[${prefix}][$(date "+%F %T")] +++ Start processing sample ${prefix} with ${sw_id}" >> ${LOG_FILE}

echo "[${prefix}][$(date "+%F %T")] 1. STARTED mapping with ${prefix}, convert, sort bam" >> ${LOG_FILE}

idx="${sw_id}_${ref_id}"
if [ "${sw_id}" = "bwa"];then
	${bwa}/bwa mem -M -t ${threads} -R "@RG\tID:$sample_id\tSM:$sample_id\tPL:Illumina" ${!idx} ${fastq1} ${fastq2} | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${threads} -m ${maxmem} -o ${work_dir}/${prefix}/${prefix}.bam -T ${tmp_dir}/${prefix}
elif [ "${sw_id}" = "gsnap" ];then
	${gsnap}/gsnap --gunzip -D ${!idx} -d ${ref_id} -t ${threads} -A sam --npath=1 -B 5 --read-group-id ${sample_id} --read-group-name ${sample_id} --read-group-platform Illumina ${fastq1} ${fastq2} | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${threads} -m ${maxmem} -o ${work_dir}/${prefix}/${prefix}.bam -T ${tmp_dir}/${prefix}
elif [ "${sw_id}" = "bowtie2" ];then
	${bowtie2}/bowtie2 -p ${threads} -x ${!idx} -1 ${fastq1} -2 ${fastq2} --rg-id "${sample_id}" --rg "SM:${sample_id}" --rg "PL:Illumina" | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${threads} -m ${maxmem} -o ${work_dir}/${prefix}/${prefix}.bam -T ${tmp_dir}/${prefix}
elif [ "${sw_id}" = "soap2" ];then
	${soap2}/soap -p ${threads} -a ${fastq1} -b ${fastq2} -D ${!idx} -o ${work_dir}/${prefix}/${prefix}_pe_output -2 ${work_dir}/${prefix}/${prefix}_se_output -u ${work_dir}/${prefix}/${prefix}_unmapped_output
	perl ${soap2sam} -p -s ${work_dir}/${prefix}/${prefix}_se_output -o ${work_dir}/${prefix}/${prefix}.sam -x ${sample_id} ${work_dir}/${prefix}/${prefix}_pe_output
	${samtools}/samtools view -bS ${work_dir}/${prefix}/${prefix}.sam | ${samtools}/samtools sort -@ ${threads} -m ${maxmem} -o ${work_dir}/${prefix}/${prefix}.bam -T ${tmp_dir}/${prefix}
elif [ "${sw_id}" = "isaac" ];then
	isaac_fastq=${work_dir}/${prefix}/isaac_fastq
	mkdir -p ${isaac_fastq}
	zcat ${fastq1} | perl -pe 's/^\+ERR[0-9]+.[0-9]+ [0-9]+ length=101/\+/' | gzip - > ${isaac_fastq}/lane1_read1.fastq.gz
	zcat ${fastq2} | perl -pe 's/^\+ERR[0-9]+.[0-9]+ [0-9]+ length=101/\+/' | gzip - > ${isaac_fastq}/lane1_read2.fastq.gz
	${isaac}/isaac-align --temp-directory ${temp_dir}/${prefix} -r ${!idx} -b ${isaac_fastq} -m 50 -j 6 --base-calls-format fastq-gz  --bam-header-tag "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:Illumina" -o ${work_dir}/${prefix} --verbosity 3 --input-parallel-load ${threads} --temp-parallel-load ${threads}
	ln -s ${work_dir}/${prefix}/Projects/default/default/sorted.bam ${work_dir}/${prefix}/${prefix}.bam
	ln -s ${work_dir}/${prefix}/Projects/default/default/sorted.bam.bai ${work_dir}/${prefix}/${prefix}.bam.bai
elif [ "${sw_id}" = "stampy" ];then
	python ${stampy}/stampy.py -g ${!idx} -h ${!idx} -t${threads} --readgroup=ID:${sample_id},SM:${sample_id},PL:Illumina --bamkeepgoodreads -M ${bwa_bam} | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${threads} -m ${maxmem} -o ${work_dir}/${prefix}/${prefix}.bam -T ${tmp_dir}/${prefix}
elif [ "${sw_id}" = "novoalign" ];then
	${novoalign}/novoalign -c ${threads} --mmapoff -t 20,3 --softclip 20 -d ${!idx} -f ${fastq1} ${fastq2} -i 350 50 -k -o SAM "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:Illumina" | ${samtools}/samtools view -Sb - | ${samtools}/samtools sort -@ ${threads} -m 4G -o ${work_dir}/${prefix}/${prefix}.bam -T ${tmp_dir}/${prefix}
fi

echo "[${prefix}][$(date "+%F %T")] 1. DONE mapping with ${prefix}, convert, sort bam" >> ${LOG_FILE}

echo "[${prefix}][$(date "+%F %T")] 2. STARTED mark duplicate/fixmate" >> ${LOG_FILE}

java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${picard}/MarkDuplicates.jar I=${work_dir}/${prefix}/${prefix}.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT AS=true OUTPUT=${work_dir}/${prefix}/${prefix}_dedup.bam METRICS_FILE=${work_dir}/${prefix}/${prefix}_dedup.log CREATE_INDEX=true
java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${picard}/FixMateInformation.jar I=${work_dir}/${prefix}/${prefix}_dedup.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true O=${work_dir}/${prefix}/${prefix}_dedup_fixmate.bam

echo "[${prefix}][$(date "+%F %T")] 2. DONE mark duplicate/fixmate" >> ${LOG_FILE}

if [ "${recal_realign_on}" = "yes" ]||[ "${recal_realign_on}" = "both" ];then
	echo "[${prefix}][$(date "+%F %T")] 3. STARTED realign/recalibration" >> ${LOG_FILE}

	java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${gatk3}/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt ${threads} --known ${dbsnp_b37} -R ${b37_fasta} --filter_reads_with_N_cigar -I ${work_dir}/${prefix}/${prefix}_dedup_fixmate.bam -o ${work_dir}/${prefix}/${prefix}.intervals
	java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${gatk3}/GenomeAnalysisTK.jar -T IndelRealigner -R ${b37_fasta} --filter_reads_with_N_cigar -I ${work_dir}/${prefix}/${prefix}_dedup_fixmate.bam -targetIntervals ${work_dir}/${prefix}/${prefix}.intervals --consensusDeterminationModel KNOWNS_ONLY -known ${b37_1000g} -LOD 0.4 -o ${work_dir}/${prefix}/${prefix}_realign.bam

	java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${gatk3}/GenomeAnalysisTK.jar -T BaseRecalibrator --filter_reads_with_N_cigar -I ${work_dir}/${prefix}/${prefix}_realign.bam -R ${b37_fasta} -knownSites ${dbsnp_b37} -o ${work_dir}/${prefix}/${prefix}.recal.table
	java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${gatk3}/GenomeAnalysisTK.jar -T PrintReads --filter_reads_with_N_cigar -I ${work_dir}/${prefix}/${prefix}_realign.bam -R ${b37_fasta} -BQSR ${work_dir}/${prefix}/${prefix}.recal.table -o ${work_dir}/${prefix}/${prefix}_recal_realign.bam
	echo "[${prefix}][$(date "+%F %T")] 3. DONE realign/recalibration" >> ${LOG_FILE}
fi

echo "[${prefix}][$(date "+%F %T")] 4. STARTED left-align indels" >> ${LOG_FILE}

if [ "${recal_realign_on}" = "no" ]||[ "${recal_realign_on}" = "both" ];then
	java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${gatk3}/GenomeAnalysisTK.jar -R ${b37_fasta} -T LeftAlignIndels -I ${work_dir}/${prefix}/${prefix}_dedup_fixmate.bam -o ${work_dir}/${prefix}/${prefix}_org_leftalignindels.bam
elif [ "${recal_realign_on}" = "yes" ]||[ "${recal_realign_on}" = "both" ];then
	java -Xmx${maxmem} -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=${tmp_dir}/${prefix} -jar ${gatk3}/GenomeAnalysisTK.jar -R ${b37_fasta} -T LeftAlignIndels -I ${work_dir}/${prefix}/${prefix}_recal_realign.bam -o ${work_dir}/${prefix}/${prefix}_recal_leftalignindels.bam
fi

echo "[${preifx}][$(date "+%F %T")] 4. DONE left-align indels" >> ${LOG_FILE}
echo "[${prefix}][$(date "+%F %T")] +++ End processing sample ${prefix} with ${sw_id}" >> ${LOG_FILE}
