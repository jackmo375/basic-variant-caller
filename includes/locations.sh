
## project path variables
pro_dir=/home/jack/Local/GeneMap/basic-variant-caller

src_dir=${pro_dir}/source
tls_dir=${pro_dir}/tools
dat_dir=${pro_dir}/data
sam_dir=${dat_dir}/sam
bam_dir=${dat_dir}/bam
vcf_dir=${dat_dir}/vcf
rds_dir=${dat_dir}/reads
ref_dir=${dat_dir}/references
fqc_dir=${dat_dir}/fastqc
tmp_dir=${dat_dir}/temp

pip_dir=${pro_dir}/pipelines

## alignment & variant calling tools
bowtie2=${tls_dir}/bowtie2-2.4.1-linux-x86_64
samtools=${tls_dir}/samtools-1.10/samtools
bcftools=${tls_dir}/bcftools-1.10.2/bcftools
bwa=${tls_dir}/bwa-0.7.17/bwa
gsnap=${tls_dir}/gmap-2019-09-12/build/bin
stampy=${tls_dir}/stampy-1.0.32
isaac=${tls_dir}/Isaac4/build/bin
gatk=${tls_dir}/gatk-4.1.7.0/gatk
picard=${tls_dir}/picard-2.23.1/picard.jar
simulate=${tls_dir}/simulate-0.1/simulate.pl
fastqc=${tls_dir}/fastqc/fastqc
trimmomatic=${tls_dir}/trimmomatic-0.39/trimmomatic-0.39.jar

jaccard=${tls_dir}/jaccard/jaccard.py
