#!/usr/bin/env python3

import argparse
from sys import exit

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Checks that the correct reference genome is used')
	parser.add_argument(
		'-f', '--fai', type=str, required=True,
		help='File with the name and length of the reference chromosomes')
	parser.add_argument(
		'-b', '--bai', type=str, required=True,
		help='File with the name of the chromosomes in the BAM')
	parser.add_argument(
		'-l', '--len', type=str, required=True,
		help='File with the length of the chromosomes in the BAM')
	parser.add_argument(
		'-p', '--pos', type=str, required=False, default=None,
		help='File with the chromosomes in the input.pos file')	
	return parser.parse_args()

def comparison(reference_dict, bam_dict, pos_file):
	## Compare both dicts
	for chrom in bam_dict:
		ref, bam = reference_dict.get(chrom), bam_dict.get(chrom)
		if ref == bam:
			continue
		else:
			return('False')

	## Check chrom from positions file if exists
	if pos_file == None:
		return('True')
	else:
		with open(pos_file) as positions:
			for line in positions:
				pos = line.strip()
				if pos in reference_dict and pos in bam_dict:
					continue
				else:
					return('False')
			return('True')

## Arguments
args = parse_args()

reference_dict, bam_dict = {}, {}
## Save reference data in a dict
with open(args.fai) as fai:
	for line in fai:
		col = line.strip().split()
		chrom, length = str(col[0]), int(col[1])
		reference_dict[chrom]=length

## Save bam data in a dict
with open(args.bai) as bai, open(args.len) as len_file:
	for line in zip(bai, len_file):
		chrom, length = str(line[0].strip()), int(line[1].strip())
		bam_dict[chrom]=length

print(comparison(reference_dict, bam_dict, args.pos))






