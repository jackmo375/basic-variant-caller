## output formatting
red='\033[0;31m'
green='\033[0;32m'
nc='\033[0m' # No Color

# generic useful functions
custom_call()
{
	# call tasks above with colored output
	# and terminate the workflow on errors

	printf "${green}$2${nc}"; echo
 
	$1 \
		|| { printf "${red}...failed!${nc}"; echo; exit 1; } \
		&& { printf "${green}...done."${nc}; echo;}
}

value_from_json() {
	file=$1		# path of input json file
	key=$2		# data field/key of {key, value} in the file
	
	# echo the value, removing any quotation marks if present
	echo $(sed -e 's/^"//' -e 's/"$//' <<<"$(cat $file | jq $key)")
}
