#!/usr/bin/env python3

import argparse
from sys import stdin

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Parser the indel distribution obtained with samtools stats to get the most common deletion and the most common insertion and their lenghts')
	parser.add_argument(
	'-i', '--input', type=str, required=True, default='-',
	help='Indel distribution obtained from the result of samtools stats with the command grep ^ID |cut -f 2- ')
	return parser.parse_args()

def main_function(file):
	for line in file:
		column = line.strip().split()
		lenght, n_in, n_del = int(column[0]), int(column[1]), int(column[2])
		if n_in == n_del:
			continue
		elif n_in > n_del:
			if n_in in in_dict:
				in_dict[n_in].append(lenght)
			else:
				in_dict[n_in] = [lenght]
			continue
		else:
			if n_del in del_dict:
				del_dict[n_del].append(lenght)
			else:
				del_dict[n_del] = [lenght]
			continue

	## Take the most common insertion and the biggest one
	try:
		max_in = max(in_dict)
		max_in_len = max(in_dict.get(max_in))
	except ValueError:
		max_in, max_in_len = 0, 0

	## Take the most common deletion and the biggest one
	try:
		max_del = max(del_dict)
		max_del_len = max(del_dict.get(max_del))
	except ValueError:
		max_del, max_del_len = 0, 0

	## Output
	print(max_in, max_in_len, max_del, max_del_len, sep='\t')

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Dictionaries
	in_dict = {}
	del_dict = {}

	## Take the input
	if args.input == '-':
		with stdin as distribution:
			main_function(distribution)
	else:
		with open(args.input) as distribution:
			main_function(distribution)
