#!/usr/bin/env python3

import argparse
from sys import stdin, argv
from re import findall, match
from statistics import mean
from math import ceil

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Filter mutations extracted by samtools mpileup using a tumor and normal bam. The output are three files: a VCF with the selected mutations and both the tumor and the normal pileup for these mutations.')
	parser.add_argument(
	'-i', '--input', type=str, required=True, metavar='MQ',
	help='Pileup file for tumor and normal in these order.')
	parser.add_argument(
	'-o', '--output_name', type=str, required=True, metavar='OUTPUT ID',
	help='Sample ID to appear in VCF')
	parser.add_argument(
	'-O', '--output_dir', type=str, required=True, default='.', metavar='OUTPUT DIR',
	help='Directory where safe the files')
	parser.add_argument(
	'-n', '--ND_name', type=str, required=True, metavar='NORMAL ID',
	help='Name for normal sample')
	parser.add_argument(
	'-t', '--TD_name', type=str, required=True, metavar='TUMOR ID',
	help='Name for tumor sample')
	parser.add_argument(
	'-s', '--sex', type=str, required=False, default='F', choices=['F','M'], metavar='SEX',
	help='Sex of the sample')
	parser.add_argument(
	'-l', '--reads_length', type=int, required=False, default=151, metavar='READS LENGTH',
	help='Length of the reads')
	parser.add_argument(
	'-Tcov', '--tumor_coverage', type=int, required=False, default=10, metavar='TUMOR COVERAGE',
	help='Minimum coverage for tumor to consider the position (default: %(default)s).')
	parser.add_argument(
	'-Ncov', '--normal_coverage', type=int, required=False, default=10, metavar='NORMAL COVERAGE',
	help='Minimum coverage for normal to consider the position (default: %(default)s).')
	parser.add_argument(
	'-Tmut', '--tumor_mut_reads', type=int, required=False, default=3, metavar='TUMOR MUT READS',
	help='Minimum number of mutated reads in tumor to select the position as mutated (default: %(default)s).')
	parser.add_argument(
	'-Nmut', '--normal_mut_reads', type=int, required=False, default=3, metavar='NORMAL MUT READS',
	help='Maximum number of mutated reads in normal to not filter the position (default: %(default)s).')
	parser.add_argument(
	'-c', '--contamination', type=float, required=False, default=0.05, metavar='CONTAMINATION',
	help='Percentage of tumor contamination in normal sample (default: %(default)s).')
	return parser.parse_args()

def correct_variants_list(variants):
	variants, indel_list = list(variants), []
	while '+' in variants or '-' in variants or '^' in variants or '$' in variants: #To convert the indels in . to count them only once
		index_list = (j for j, k in enumerate(variants) if k == '+' or k == '-' or k == '^' or k == '$') #j is the index and k the element in variants, so if the element is + or -, it returns its index
		n = next(index_list) #To focus on the first element
		if variants[n] == '^': #^ represents the start of the read and the following character shows the mapQ of the read, thus we delete one of the symbols to take into account only once that mapQ value
			variants[n:n+2] = ''
			continue
		elif variants[n] == '$': #$ represents the end of a read segment, but the read mapQ value doesn't appear, so we have to remove the symbol
			variants[n:n+1] = ''
			continue
		else: #indel
			try:
				size = int(str(variants[n+1])+str(variants[n+2]))+2 #In case the indel had more than 9 changes
			except ValueError:
				size = int(variants[n+1])+1 #The size of the indel, which is the element following the - or + symbol
			indel_list.append(''.join(variants[n:n+size+1]))
			variants[n-1:n+size+1] = ';' #To convert the indel pattern into a . (e.g. .+4ACGT > ;), thus there is only one symbol for each mapQ value. Note that before the indel pattern there is a . or , which is the symbol that has the mapQ value, so we also have to replace it
			continue
	variants_list = findall(r'(?<!\^)[ACGT]', ''.join(variants)) #To select the substitutions
	single_variants = set(variants_list) #Variants without duplicates
	return(single_variants, variants_list, variants)

def most_common_variant(single_variants_list, full_variants_list):
	most_count = 0
	for variant in single_variants_list:
		count = full_variants_list.count(variant)
		if count > most_count:
			most_count = count
			alt = variant
		else:
			continue
	return(alt, int(most_count))

def get_mut_position(alt, variants, positions):
	positions_list = []
	for base in zip(variants, positions):
		if base[0] == alt: #To extract only those reads that contain the mutation we are looking for
			positions_list.append(int(base[1]))
		else:
			continue
	distance = max(positions_list)-min(positions_list)
	mean_distance = mean(positions_list)/args.reads_length
	return(distance, mean_distance)

def read_names(alt, variants, read_names):
	names_list = []
	for base in zip(variants, read_names):
		if base[0] == alt: #To extract only those reads that contain the mutation we are looking for
			names_list.append(base[1])
		else:
			continue
	return(names_list)

def decode_ASCII(code):
	decode_value = []
	for i in code:
		decode_value.append(ord(i)-33) # ASCII code -33 (samtools)
	return(float(mean(decode_value)))

def main_function(pileup):
	for line in pileup:
		column = line.strip().split()
		chrom, pos, ref, cov_tum, variants_tum, baseq_tum, positions_tum, read_names_tum, mapping_tum, cov_norm, variants_norm, baseq_norm, positions_norm, mapping_norm = str(column[0]), column[1], column[2], int(column[3]), column[4].upper(), column[5], column[6].split(','), column[7].split(','), column[8], int(column[9]), column[10].upper(), column[11], column[12].split(','), column[13]

		## Set the coverage threshold depending on the sample's sex
		if args.sex == 'M' and (chrom == 'X' or chrom == 'Y'):
			tumor_coverage = float(args.tumor_coverage/2)
			normal_coverage = float(args.normal_coverage/2)
		elif args.sex == 'F' and chrom == 'Y':
			continue
		else:
			tumor_coverage = args.tumor_coverage
			normal_coverage = args.normal_coverage
		
		## Filter by coverage
		if cov_tum >= tumor_coverage and cov_norm >= normal_coverage:
			single_variants_tum, variants_list_tum, new_variants_tum = correct_variants_list(variants_tum)
			single_variants_norm, variants_list_norm, new_variants_norm = correct_variants_list(variants_norm)

			## Check the number of variants for both tumor and normal
			if len(single_variants_tum) == 0: # If there isn't any mutation in the tumor
				continue
			else:
				alt_tum, alt_count_tum = most_common_variant(single_variants_tum, variants_list_tum)
				if len(single_variants_norm) == 0: # If there isn't any mutation in the normal
					alt_norm, alt_count_norm = '.', 0
				else:
					alt_count_norm = variants_list_norm.count(alt_tum)
					if alt_count_norm > 0:
						alt_norm = alt_tum
					else:
						alt_norm = '.'

			## Check in which position of the reads the mutation appears
			distance_position, mean_position = get_mut_position(alt_tum, new_variants_tum, positions_tum)

			## Remove those mutations that appear less than N times and mutations that also appear in the normal
			# Calculate tumor's MAF
			MAF_tum = alt_count_tum/cov_tum
			# Apply contamination to the max number of mut reads in normal
			max_normal_mut_reads = ceil(args.normal_mut_reads + (cov_norm*args.contamination*MAF_tum))
			# Filter
			if alt_count_tum < args.tumor_mut_reads or (alt_tum == alt_norm and alt_count_norm > max_normal_mut_reads):
				continue
			else:
				if len(alt_tum) == 1:
					pass
				else:
					nucleotides = ''.join(findall(r'[A-Z]', alt_tum))
					indel_deletion = match('-', alt_tum)
					if indel_deletion is not None:
						alt_tum = ref
						ref += nucleotides
					else:
						alt_tum = ref + nucleotides

				## Save the name of the mutated reads in tumor
				read_names_list = read_names(alt_tum, new_variants_tum, read_names_tum)

				## Calculate the MAF, BQ and MQ for the position of the mutation
				MAF = alt_count_tum/cov_tum
				BQN, BQT = decode_ASCII(baseq_norm), decode_ASCII(baseq_tum)
				MQN, MQT = decode_ASCII(mapping_norm), decode_ASCII(mapping_tum)

				## Write the output
				print(chrom, pos, args.output_name, ref.upper(), alt_tum.upper(), '.', '.', 'VAF=%.3f' % float(MAF), 'DP:BQ:MQ:AD', '%i:%.1f:%.1f:%i' % (cov_norm, BQN, MQN, alt_count_norm), '%i:%.1f:%.1f:%i' % (cov_tum, BQT, MQT, alt_count_tum), sep='\t', end='\n')
				ND_pileup.write(chrom + "\t" + pos + "\t" + ref + "\t" + str(cov_norm) + "\t" + variants_norm + "\t" + baseq_norm + "\t" + mapping_norm + "\t" + str(','.join(positions_norm)) + "\n")
				TD_pileup.write(chrom + "\t" + pos + "\t" + ref + "\t" + str(cov_tum) + "\t" + variants_tum + "\t" + baseq_tum + "\t" + mapping_tum + "\t" + str(','.join(positions_tum)) + "\n")
				mutations_interval.write(chrom + '_' + pos + "\t" +  str(distance_position) + "\t" + str(mean_position) + "\n")
				read_names_file.write(chrom + '_' + pos + "\t" + ','.join(read_names_list) + "\n")
		else:
			continue

	## Close output files
	ND_pileup.close()
	TD_pileup.close()


if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Define the mini.pileup outputs for the ML pipeline
	ND_pileup = open(args.output_dir+'/'+args.ND_name+'.mini.pileup', 'w+')
	TD_pileup = open(args.output_dir+'/'+args.TD_name+'.mini.pileup', 'w+')
	mutations_interval = open(args.output_dir+'/'+args.output_name+'.mutations.interval', 'w+')
	read_names_file = open(args.output_dir+'/'+args.TD_name+'.read.names', 'w+')

	## VCF header
	arg_list = []
	for arg in sorted(vars(args)):
		arg_list.append('--'+arg)
		arg_list.append(str(getattr(args, arg)))
	header = ['#CHROM','POS','ID','REF','ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.ND_name, args.TD_name]
	print('##fileformat=VCFv4.2', '##INFO=<ID=VAF,Number=1,Type=Float,Description="Allele frequency of the alternative allele">',  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth coverage">',	'##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average base quality">', '##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Average mapping quality">', '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Number of alternative reads">', '##Command=python3 %s %s' % (argv[0], ' '.join(arg_list)), '\t'.join(header), sep='\n')

	mutations_interval.write('ID' + "\t" + 'INTERVAL_SIZE' + "\t" + 'MEAN_POSITION' + "\n")

	## Take the input
	if args.input == '-':
		with stdin as pileup:
			main_function(pileup)
	else:
		with open(args.input) as pileup:
			main_function(pileup)
