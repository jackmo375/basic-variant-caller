#
#	GENERATE JSON - python tool
#
#	* generate an input json file series
#	* to input to pipelines (for testing
#	* the variation of pipeline output
#	* with input parameters). 
#
#	Jack Morrice
#
##########################################

import sys, argparse
import json


def main():

	args=get_shell_arguments()

	with open(args.input_json) as json_file:

		# read input json file
		inputs = json.load(json_file)

		# modify the chosen value
		inputs[args.key] = args.new_value

		# write out to file
		output=open(args.output_json, 'w')
		output.write(json.dumps(inputs, indent=2))
		output.close()


#
#  functions
#
def get_shell_arguments():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--input_json',
		type=str,
		help='input json file to modify')

	parser.add_argument(
		'--key',
		type=str,
		help='key of value to modify')

	parser.add_argument(
		'--new_value',
		type=str,
		help='new value to replace the input one at <key> location')

	parser.add_argument(
		'--output_json',
		type=str,
		help='output json file name to save to')

	return parser.parse_args()


if __name__ == '__main__':
	main()