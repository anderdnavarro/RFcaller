#!/usr/bin/env python3

import argparse
from pysam import faidx
from regex import findall #We use this module instead of the python's re module because with this we can search overlapping matches with the findall option

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script counts how many dimers of each type has a specific sequence')
	parser.add_argument(
	'-p', '--position', type=str, required=True, metavar='POSITIONS_BED',
	help='Bed file containing the positions to analyze')
	parser.add_argument(
	'-w', '--window', type=int, required=False, default=15, metavar='WINDOW',
	help='Window size (in base pairs) to expand the position (default: %(default)s)')
	parser.add_argument(
	'-r', '--ref_path', type=str, required=True, metavar='REFERENCE_GENOME',
	help='Path to the reference genome')
	return parser.parse_args()

def main_function(file):
	for line in file:
		column = line.strip().split()
		chrom, pos = str(column[0]), int(column[1])
		print(chrom, pos, sep='\t', end='\t')

		start, end = str(pos - args.window), str(pos + args.window) #To extend the region
		region = chrom + ':' + start + '-' + end
		sequence = faidx(args.ref_path, region).splitlines() #To search the sequence in the reference genome

		values_list = []
		for dimer in dimers_list:
			dimer_search = findall(dimer, sequence[1], overlapped=True)
			print(len(dimer_search), end='\t')
			values_list.append(len(dimer_search))
		
		## Calculate %GC of the region
		GC = (len(findall('G', sequence[1]))+len(findall('C', sequence[1])))/(int(end)-int(start))
		## Print results
		print(max(values_list), GC, sep='\t')


if __name__ == '__main__':
	## Arguments
	args = parse_args()

	# Dimer list
	dimers_list = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']

	# Header
	print('#Chr', 'Pos', '\t'.join(dimers_list), 'Max_Value', 'GC_percentaje', sep='\t', end='\n')

	with open(args.position) as pos_file:
		main_function(pos_file)
