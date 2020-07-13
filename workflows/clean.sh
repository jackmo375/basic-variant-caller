
source ../configs/parameters.cfg

echo "cleaning bam folder" && rm -r ${bam_dir}/* || exit 1
echo "cleaning sam folder" && rm -r ${sam_dir}/* || exit 1
echo "cleaning vcf folder" && rm -r ${vcf_dir}/* || exit 1
echo "cleaning index folder" && rm -r ${dat_dir}/indices/* || exit 1
echo "cleaning references folder" && rm -r ${dat_dir}/references/*.fai || exit 1
echo "cleaning temp folder" && rm -r ${tmp_dir}/* || exit 1
