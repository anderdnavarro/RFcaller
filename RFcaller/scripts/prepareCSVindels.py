#!/usr/bin/env python3

import argparse
from re import findall

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script is used to make a CSV file that can be used by the machine learning algorithm to detect mutations.')
	parser.add_argument(
	'-n', '--name', type=str, required=True, metavar='NAME',
	help='output name')
	parser.add_argument(
	'--ND_pileup', type=str, required=True, metavar='NORMAL PILEUP',
	help='mini.pileup for normal sample')
	parser.add_argument(
	'--TD_pileup', type=str, required=True, metavar='TUMOR PILEUP',
	help='mini.pileup for tumor sample')
	parser.add_argument(
	'--ND_mapq', type=str, required=True, metavar='NORMAL MAPQ',
	help='pileup with mapping quiality for normal sample')
	parser.add_argument(
	'--TD_mapq', type=str, required=True, metavar='TUMOR MAPQ',
	help='pileup with mapping quiality for tumor sample')
	parser.add_argument(
	'-s', '--stats', type=str, required=True, metavar='STATS',
	help='file with BAM statistics')
	parser.add_argument(
	'-c', '--cigar', type=str, required=True, metavar='CIGAR',
	help='file with cigar information')
	parser.add_argument(
	'-d', '--sequence', type=str, required=True, metavar='SEQUENCE',
	help='file with information about the context (sequence)')
	parser.add_argument(
	'-i', '--interval', type=str, required=True, metavar='INTERVAL',
	help='file with information about the position of mutations along the reads (interval)')
	parser.add_argument(
	'-dis', '--distribution', type=str, required=True, metavar='DISTRIBUTION',
	help='file with information about the distribution of the indels')
	parser.add_argument(
	'-f', '--features', type=str, required=True, metavar='FEATURES',
	help='file with information about some features of the indels')
	return parser.parse_args()

## To extract each field for each file
def start(file, file_format):
	column = file.strip().split()
	if file_format == 'pileup':
		chrom, pos, coverage, bases = str(column[0]), str(column[1]), int(column[3]), column[4].upper()
		return(chrom, pos, coverage, bases)
	elif file_format == 'mapq':
		chrom, pos, coverage, bases, mapQ = str(column[0]), str(column[1]), int(column[3]), column[4].upper(), column[6]
		return(chrom, pos, coverage, bases, mapQ)
	elif file_format == 'sequence':
		return(str(column[-2]), str(column[-1]))
	elif file_format == 'interval':
		return(str(column[1]), str(column[2]))
	elif file_format == 'features':
		return(str(column[1]))
	elif file_format == 'distribution':
		return(column[1:])

def correct_variants_list(variants):
	variants, indel_list = list(variants), []
	while '+' in variants or '-' in variants or '^' in variants or '$' in variants: #To convert the indels in ; to count them only once
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
			variants[n-1:n+size+1] = ';' #To convert the indel pattern into a . (e.g. .+4ACGT > ;), thus there is only one symbol for each read. Note that before the indel pattern there is a . or , which is the symbol that has the mapQ value, so we also have to replace it
			continue
	## Get the most frequent variant
	if len(indel_list) == 0:
		alt_indel, indel_count = '', 0
	else:
		alt_indel, indel_count = most_common_indel(set(indel_list), indel_list)
	return(variants, alt_indel, indel_count)

def most_common_indel(single_indel_list, full_indel_list):
	most_count = 0
	for indel in single_indel_list:
		count = full_indel_list.count(indel)
		if count >= most_count:
			most_count = count
			alt = indel
		else:
			continue
	return(alt, most_count)

## To calculate the error_ratio of the region
def extract_error_ratio(file, normal_mut, tumor_mut):
	column = file.strip().split()
	normal_base, normal_mismatch, tumor_base, tumor_mismatch = int(column[1]), int(column[2]), int(column[3]), int(column[4])
	normal_error_ratio = (normal_mismatch-normal_mut)/(normal_base-normal_mut)
	tumor_error_ratio = (tumor_mismatch-tumor_mut)/(tumor_base-tumor_mut)
	return(normal_error_ratio, tumor_error_ratio)

## To process the cigar field
def extract_cigar(file):
	column = file.strip().split()
	bad_normal_cigar = findall(r'([HS])', column[1])
	bad_tumor_cigar = findall(r'([HS])', column[2])
	return(len(bad_normal_cigar), len(bad_tumor_cigar))

def prepare_mapQ(variants, mapQ_column, alt):
	base_mapQ_list = list(mapQ_column) #Extract information from the mq
	if len(variants) == len(base_mapQ_list): #The number of symbols in the variants list has to be the same as the number of mapQ values
		mut_mapQ, normal_mapQ = [],[]
		for base in zip(variants, base_mapQ_list): #Now we associate each base with each mapQ
			if base[0] == ';': #To extract only those read names that contain the mutation we are looking for; #MUT base
				map_value = ord(base[1]) - 33 #Samtools says that the mapQ is the ASCII number - 33
				mut_mapQ.append(map_value)
			else: #WT base
				map_value = ord(base[1]) - 33 #Samtools says that the mapQ is the ASCII number - 33
				normal_mapQ.append(map_value)
		norm_reads_number, mut_reads_number = len(normal_mapQ), len(mut_mapQ)
		try:
			normal_reads = sum(normal_mapQ)/norm_reads_number
		except ZeroDivisionError:
			normal_reads = 0
		try:
			mut_reads = sum(mut_mapQ)/mut_reads_number
		except ZeroDivisionError:
			mut_reads = 0
		return(norm_reads_number, mut_reads_number, normal_reads, mut_reads)

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Header
	print('#ID', 'Q30N_cov', 'Q30T_cov', 'Q30N_mut_reads', 'Q30T_mut_reads', 'N_mut_reads', 'T_mut_reads', 'N_normal_mapQ', 'T_normal_mapQ', 'N_mut_mapQ', 'T_mut_mapQ', 'Normal_Error_Ratio', 'Tumor_Error_Ratio', 'Normal_cigar', 'Tumor_cigar', 'Dimers', 'GC_percentaje', 'Interval_size', 'Mean_position', 'N_repeat_indel', 'Normal_insertion_count', 'Normal_insertion_lenght', 'Normal_deletion_count', 'Normal_deletion_lenght', 'Tumor_insertion_count', 'Tumor_insertion_lenght', 'Tumor_deletion_count', 'Tumor_deletion_lenght', sep=',', end='\n')
	
	## Open files
	Q30N, Q30T, TN, TT, stats, cigar, sequence, interval, features, distribution = open(args.ND_pileup), open(args.TD_pileup), open(args.ND_mapq), open(args.TD_mapq), open(args.stats), open(args.cigar), open(args.sequence), open(args.interval), open(args.features), open(args.distribution)

	## Skip headers
	next(stats)
	next(cigar)
	next(sequence)
	next(interval)
	next(features)
	next(distribution)

	for line in zip(Q30N, Q30T, TN, TT, stats, cigar, sequence, interval, features, distribution): #We read the ten files at the same time, because they have the same mutation for each line
		try:
			chrom, pos, Q30N_cov, Q30N_bases = start(line[0], 'pileup')
			chrom, pos, Q30T_cov, Q30T_bases = start(line[1], 'pileup')
			chrom, pos, TN_cov, TN_bases, TN_mapQ = start(line[2], 'mapq')
			chrom, pos, TT_cov, TT_bases, TT_mapQ = start(line[3], 'mapq')

			## Correct the variants field and take the alternative base for each MQ
			Q30N_bases_correct, alt_Q30N, mut_reads_Q30N = correct_variants_list(Q30N_bases)
			Q30T_bases_correct, alt_Q30T, mut_reads_Q30T = correct_variants_list(Q30T_bases)
			TN_bases_correct, alt_TN, mut_reads_TN = correct_variants_list(TN_bases)
			TT_bases_correct, alt_TT, mut_reads_TT = correct_variants_list(TT_bases)
			## Filter positions with many mut reads in normal with low quality
			if (mut_reads_TN-mut_reads_Q30N > 0.05*TN_cov) and (mut_reads_TN/TN_cov > 0.1):
				continue
			else:
				norm_reads_number_TN, mut_reads_number_TN, normal_mapQ_TN, mut_mapQ_TN = prepare_mapQ(TN_bases_correct, TN_mapQ, alt_TT)
				norm_reads_number_TT, mut_reads_number_TT, normal_mapQ_TT, mut_mapQ_TT = prepare_mapQ(TT_bases_correct, TT_mapQ, alt_TT)
				normal_error_ratio, tumor_error_ratio = extract_error_ratio(line[4], mut_reads_number_TN, mut_reads_number_TT)
				normal_cigar, tumor_cigar = extract_cigar(line[5])
				max_dimer_value, GC_percentaje = start(line[6], 'sequence')
				distance_position, mean_position = start(line[7], 'interval')
				features_values = start(line[8], 'features')
				distribution_values = start(line[9], 'distribution')
				print('%s_%s_%s' % (chrom, pos, args.name), Q30N_cov, Q30T_cov, mut_reads_Q30N, mut_reads_Q30T, mut_reads_TN, mut_reads_TT, normal_mapQ_TN, normal_mapQ_TT, mut_mapQ_TN, mut_mapQ_TT, normal_error_ratio, tumor_error_ratio, normal_cigar, tumor_cigar, max_dimer_value, GC_percentaje, distance_position, mean_position, features_values, ','.join(distribution_values), sep=',', end='\n')
		except ZeroDivisionError:
			continue

	Q30N.close()
	Q30T.close()
	TN.close()
	TT.close()
	stats.close()
	cigar.close()
	sequence.close()
	interval.close()
	features.close()
	distribution.close()
