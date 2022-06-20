#!/usr/bin/env python3

import argparse
from sys import stdin
from pysam import faidx

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Script to gather mutations that are together in the same strand')
	parser.add_argument(
		'-i', '--input', type=str, required=True, metavar='VCF',
		help='Sorted VCF')
	parser.add_argument(
		'-r', '--read_names', type=str, required=True, metavar='READ_NAMES',
		help='File with mutated read names in tumor')
	parser.add_argument(
		'-f', '--fasta', type=str, required=True, metavar='GENOME',
		help='Reference genome in fasta format')
	return parser.parse_args()

def check_dinucleotides(chrom, prev_pos, pos):
	prev_names, actual_names = read_names_dict[str(chrom)+'_'+str(prev_pos)], read_names_dict[str(chrom)+'_'+str(pos)]
	coincidences, min_list = len([i for i in prev_names if i in actual_names]), min([len(prev_names),len(actual_names)])
	coincidence_freq = float(coincidences/min_list)
	return(coincidence_freq)

def main_function(file):
	first = True
	for line in file:
		## Print header
		if line.startswith('#'):
			print(line.strip())
			continue
		else:
			column = line.strip().split()
			chrom, pos, ref, alt, AF = str(column[0]), int(column[1]), str(column[3]), str(column[4]), float(column[7].split(';')[0].split('=')[1])
			## Save first mutation
			if first:
				prev_line = line.strip()
				prev_chrom, prev_pos, prev_ref, prev_alt, prev_AF = chrom, pos, ref, alt, AF
				first = False
				continue
			else:
				## Check whether the next mutation is in the next position
				if prev_chrom == chrom and (pos-prev_pos) == 1:
					coincidence_freq = check_dinucleotides(chrom, prev_pos, pos)
					if coincidence_freq >= 0.2:
						prev_col = prev_line.split()
						prev_col[3], prev_col[4] = (prev_ref+ref), (prev_alt+alt)
						print('\t'.join(prev_col))
						first = True
					else:
						print(prev_line)
						prev_line = line.strip()
						prev_chrom, prev_pos, prev_ref, prev_alt, prev_AF = chrom, pos, ref, alt, AF
					continue
				elif prev_chrom == chrom and (pos-prev_pos) == 2:
					coincidence_freq = check_dinucleotides(chrom, prev_pos, pos)
					if coincidence_freq >= 0.2:
						mid_base = faidx(args.fasta, chrom+':'+str(prev_pos+1)+'-'+str(prev_pos+1)).splitlines()[1]
						prev_col = prev_line.split()
						prev_col[3], prev_col[4] = (prev_ref+str(mid_base)+ref), (prev_alt+str(mid_base)+alt)
						print('\t'.join(prev_col))
						first = True
					else:
						print(prev_line)
						prev_line = line.strip()
						prev_chrom, prev_pos, prev_ref, prev_alt, prev_AF = chrom, pos, ref, alt, AF
					continue
				else:
					print(prev_line)
					prev_line = line.strip()
					prev_chrom, prev_pos, prev_ref, prev_alt, prev_AF = chrom, pos, ref, alt, AF
					continue
	if not first:
		print(prev_line)

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Read names dictionary
	read_names_dict = {}
	with open(args.read_names) as names:
		for line in names:
			col = line.strip().split()
			read_names_dict[col[0]] = col[1].split(',')

	## Read input
	if args.input == '-':
		with stdin as file:
			main_function(file)
	else:
		with open(args.input) as file:
			main_function(file)
