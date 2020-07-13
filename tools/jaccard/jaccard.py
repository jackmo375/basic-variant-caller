#
#  *Jaccard.py*
#	compute jaccard index of vcf files
#	output from: "bcftools isec"
#
########################################

import sys, argparse

def main():
	args = getShellArguments()

	instream = open(args.input, 'r')

	lines = instream.readlines()

	# check both intersections are the same:
	if (int(lines[5].strip()) == int(lines[7].strip())):
		n_setdiff = [int(lines[1].strip()), int(lines[3].strip())]
		n_isec    = int(lines[5].strip())
	else:
		raise Exception('set intersection files are not the same')

	J = n_isec / (n_isec + sum(n_setdiff))

	print("******************")
	print("* Jaccard index: *")
	print("{0:8.3f}".format(J))
	print("******************")


#
#  functions
#
def getShellArguments():
	'''
	Get shell enviroment arguments
	'''
	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--input',
		type=str,
		help='path to bcftools isec output file summary')

	return parser.parse_args()

#
#  run main
#
if __name__ == '__main__':
	main()
