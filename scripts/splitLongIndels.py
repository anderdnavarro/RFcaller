#!/usr/bin/env python3

import argparse

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Split long and short indels into two different files')
	parser.add_argument(
		'-i', '--input', type=str, required=True, metavar='VCF',
		help='VCF with INDELs')
	parser.add_argument(
		'-s', '--size', type=int, required=True, metavar='SIZE',
		help='Size of the indel to consider it as long')
	parser.add_argument(
		'-b', '--baseName', type=str, required=True, metavar='NAME',
		help='Output basename')	
	parser.add_argument(
		'-o', '--outDir', type=str, required=True, metavar='DIR',
		help='Output directory')
	return parser.parse_args()

## Arguments
args = parse_args()

## Output files
short_indels = open(args.outDir+'/short_'+args.baseName, '+w')
long_indels = open(args.outDir+'/long_'+args.baseName, '+w')
all_indels = open(args.outDir+'/'+args.baseName, '+w')

## Main function
with open(args.input) as file:
	for line in file:
		if line.startswith('#'):
			continue
		else:
			col = line.strip().split()
			chrom, pos, ref, alt = str(col[0]), str(col[1]), str(col[3]), str(col[4])
			if len(ref) >= args.size or len(alt) >= args.size:
				long_indels.write(chrom+'\t'+pos+'\n')
			else:
				short_indels.write(chrom+'\t'+pos+'\n')
			all_indels.write(chrom+'\t'+pos+'\n')
			continue
			
