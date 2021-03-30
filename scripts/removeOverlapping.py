#!/usr/bin/env python3

import argparse
from sys import stdin

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Remove overlapping reads from a SAM')
	parser.add_argument(
	'-i', '--input', type=str, required=False, default='-', metavar='SAM',
	help='Read sorted SAM file')
	return parser.parse_args()

def main_function(file):
	for line in file:
		if line.startswith('@'): #Skip header
			print(line.strip())
			continue
		else:
			column = line.strip().split()
			read_ID = '_'.join(column[:2])
			if read_ID not in reads_list: #Print only single reads
				print(line.strip())
				if len(reads_list) == 50:
					del reads_list[0]
					reads_list.append(read_ID)
				else:
					reads_list.append(read_ID)
			else:
				continue

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	reads_list = []
	## Take the input
	if args.input == '-':
		with stdin as bam:
			main_function(bam)
	else:
		with open(args.input) as bam:
			main_function(bam)
