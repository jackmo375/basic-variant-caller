## output formatting
red='\033[0;31m'
green='\033[0;32m'
nc='\033[0m' # No Color
error="\e[38;5;1mERROR\e[0m:"

# generic useful functions

custom_call() {
	# for *top level functions*
	# call tasks above with colored output
	# and terminate the workflow on errors

	printf "${green}$2${nc}"; echo
 
	$1 \
		|| { printf "${red}...failed!${nc}"; echo; exit 1; } \
		&& { printf "${green}...done."${nc}; echo; }
}

value_from_json() {
	local file=$1		# path of input json file
	local key=$2		# data field/key of {key, value} in the file
	local dest=$3		# destination to save value to
	
	# save the value, removing any quotation marks if present
	local value=$(sed -e 's/^"//' -e 's/"$//' <<<"$(cat $file | jq $key)")

	# if the value is not null then the key exists in the json file, 
	# and we can update the $dest variable with $value
	[[ $value != "null" ]] && eval $dest=$value
}

random_id() {
	openssl rand -hex 4
}

check_int() {
	local value=$1
	local name=$2

	re=^[0-9]+$
	if ! [[ $1 =~ ^[0-9]+$ ]] ; then
   		echo "error: $2 is not a positive integer" >&2; return 1
	fi
}