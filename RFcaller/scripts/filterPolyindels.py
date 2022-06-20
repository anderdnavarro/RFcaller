#!/usr/bin/env python3

import argparse
from sys import stdin

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Script to choose indels of the type: CAAAAA > CAAAA')
	parser.add_argument(
		'-i', '--input', type=str, required=False, metavar='VCF', default = '-',
		help='VCF file from initial calling')
	return parser.parse_args()

def main_function(file):
	for line in file:
		if line.startswith('#'):
			continue
		else:
			col = line.strip().split()
			chrom, pos, ref, alt, qual = str(col[0]), int(col[1]), str(col[3]), str(col[4]), float(col[5])
			if (len(ref) > len(alt)) and (abs(len(ref)-len(alt))<=2):
				exp_ref = alt+alt[-1]
				if exp_ref == ref:
					print(chrom, pos, pos, sep='\t' , end='\n')
			elif len(alt) > len(ref) and (abs(len(ref)-len(alt))<=2):
				exp_alt = ref+ref[-1]
				if exp_alt == alt:
					print(chrom, pos, pos, sep='\t' , end='\n')
			else:
				continue

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Take the input
	if args.input == '-':
		with stdin as file:
			main_function(file)
	else:
		with open(args.input) as file:
			main_function(file)


