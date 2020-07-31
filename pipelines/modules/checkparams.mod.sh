#
#  checkparams (module)
#
#	summary:
#   + check_input_json()
#   + check_id()
#	+ check_ref()
#	+ check_bwa_idx()
#	+ check_gatk_dict()
#	+ check_samtools_fai()
#
check_input_json() {

    # check a filename was given
    [[ -s ${inputs["input_json"]} ]] || { echo 'ERROR: please specify an input json file'; return 1; }

    # check that it is a valid json file
    cat ${inputs["input_json"]} | jq -e . 1> /dev/null || { echo 'ERROR in input json file'; return 1; }
}

check_id() {
    local \
        id_type=$1 \
        id=$2 \
        status=0

    # check the id is not empty or null
    [[ -z $id ]] && { echo 'ERROR: input $id_type id is empty, please choose another'; status=1; }

    # check the id type has the relevant identifier: c_..., s_...
    if [[ $id_type == 'cohort' ]]; then
        [[ $id == "c_"* ]] || { echo "ERROR: cohort ids must start with c_... please choose another"; status=1; }
    elif [[ $id_type == "s_"* ]]; then
        [[ $id == "s_"* ]] || { echo "ERROR: sample ids must start with s_... please choose another"; status=1; }
    fi

    # check the id has no spaces or tabs

    return $status
}

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