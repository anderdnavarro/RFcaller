#!/usr/bin/env python3

import argparse
from re import findall
import numpy as np
import pysam

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script is used to make a CSV file that can be used by the machine learning algorithm to detect mutations.')
	parser.add_argument(
	'-n', '--name', type=str, required=True, metavar='NAME',
	help='output name')
	parser.add_argument(
	'-b', '--tumor_bam', type=str, required=True, metavar='TUMOR_BAM',
	help='reduced bam for tumor sample')
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
	help='file with information about the position of mutations along the reads')
	parser.add_argument(
	'-r', '--read_names', type=str, required=True, metavar='READ_NAMES',
	help='file with information about the name of the mutated reads for each position')
	parser.add_argument(
	'-Tmut', '--tumor_mut_reads', type=int, required=False, default=3, metavar='TUMOR MUT READS',
	help='Minimum number of mutated reads in tumor to select the position as mutated (default: %(default)s).')
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

def most_common_variant(single_variants_list, full_variants_list):
	most_count = 0
	for variant in single_variants_list:
		count = full_variants_list.count(variant)
		if count >= most_count:
			most_count = count
			alt = variant
		else:
			continue
	return(alt, most_count)

def correct_variants_list(variants):
	variants, alts_list = list(variants.upper()), []

	## Correct the variant field
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
			alts_list.append(''.join(variants[n:n+size+1]))
			variants[n-1:n+size+1] = ';' #To convert the indel pattern into a . (e.g. .+4ACGT > ;), thus there is only one symbol for each mapQ value. Note that before the indel pattern there is a . or , which is the symbol that has the mapQ value, so we also have to replace it
			continue

	## Get the most frequent variant
	alts_list += findall(r'(?<!\^)[ACGT]', ''.join(variants))
	single_variants = set(alts_list) #Variants without duplicates
	if len(single_variants) == 0: #This can occur if the read is WT
		alt_count, indel_alt, alt = 0, False, ''
	else:
		alt, alt_count = most_common_variant(single_variants, alts_list) #To obtain the most common variant and how many times it appears
		if len(alt) > 1:
			indel_alt = True
		else:
			indel_alt = False
	return(variants, alt, alt_count, indel_alt)

## To calculate the error_ratio of the region
def extract_error_ratio(file, normal_mut, tumor_mut):
	line = file.strip()
	column = line.split()
	normal_base, normal_mismatch, tumor_base, tumor_mismatch = int(column[1]), int(column[2]), int(column[3]), int(column[4])
	normal_error_ratio = (normal_mismatch-normal_mut)/(normal_base-normal_mut)
	tumor_error_ratio = (tumor_mismatch-tumor_mut)/(tumor_base-tumor_mut)
	return(normal_error_ratio, tumor_error_ratio)

## To process the cigar field
def extract_cigar(file):
	line = file.strip()
	column = line.split()
	bad_normal_cigar = findall(r'([HSID])', column[1])
	bad_tumor_cigar = findall(r'([HSID])', column[2])
	return (len(bad_normal_cigar), len(bad_tumor_cigar))

def prepare_mapQ(variants, mapQ_column, alt):
	base_mapQ_list = list(mapQ_column) #Extract information from the mq
	if len(variants) == len(base_mapQ_list): #The number of symbols in the variants list has to be the same as the number of mapQ values
		mut_mapQ, normal_mapQ = [],[]
		for base in zip(variants, base_mapQ_list): #Now we associate each base with each mapQ
			if base[0] == alt: #To extract only those read names that contain the mutation we are looking for; #MUT base
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

def baseq_filter(samfile, reads_file, threshold = 15):
	passed_reads = 0
	chrom_pos, reads = reads_file.strip().split()
	chrom, pos = chrom_pos.split('_')
	region = str(chrom) + ":" + str(pos) + "-" + str(pos)
	for read in samfile.fetch(region = region):
		if read.query_name in reads:
			quals = [ord(x) - 33 for x in read.qual]     
			Q1 = np.quantile(quals, .25)
			
			if Q1 > threshold:
				passed_reads += 1
		else:
			continue
	return passed_reads

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Header
	print('#ID', 'Q30N_cov', 'Q30T_cov', 'Q30N_mut_reads', 'Q30T_mut_reads', 'N_normal_mapQ', 'T_normal_mapQ', 'N_mut_mapQ', 'T_mut_mapQ', 'Normal_Error_Ratio', 'Tumor_Error_Ratio', 'Normal_cigar', 'Tumor_cigar', 'Dimers', 'GC_percentage', 'Interval_size', 'Mean_position', sep=',', end='\n')

	## Open files
	Q30N, Q30T, TN, TT, stats, cigar, sequence, interval, names, samfile = open(args.ND_pileup), open(args.TD_pileup), open(args.ND_mapq), open(args.TD_mapq), open(args.stats), open(args.cigar), open(args.sequence), open(args.interval), open(args.read_names), pysam.AlignmentFile(args.tumor_bam, "rb")

	## Skip headers
	next(stats)
	next(cigar)
	next(sequence)
	next(interval)

	for line in zip(Q30N, Q30T, TN, TT, stats, cigar, sequence, interval, names): #We read the seven files at the same time, because they have the same mutation for each line
		try:
			chrom, pos, Q30N_cov, Q30N_bases = start(line[0], 'pileup')
			chrom, pos, Q30T_cov, Q30T_bases = start(line[1], 'pileup')
			chrom, pos, TN_cov, TN_bases, TN_mapQ = start(line[2], 'mapq')
			chrom, pos, TT_cov, TT_bases, TT_mapQ = start(line[3], 'mapq')

			## Correct the variants field and take the alternative base for each MQ
			Q30N_bases_correct, alt_Q30N, mut_reads_Q30N, Q30N_indel = correct_variants_list(Q30N_bases)
			Q30T_bases_correct, alt_Q30T, mut_reads_Q30T, Q30T_indel = correct_variants_list(Q30T_bases)
			TN_bases_correct, alt_TN, mut_reads_TN, TN_indel = correct_variants_list(TN_bases)
			TT_bases_correct, alt_TT, mut_reads_TT, TT_indel = correct_variants_list(TT_bases)
			if (Q30N_indel and (mut_reads_Q30N/Q30N_cov) > 0.1) or (TN_indel and (mut_reads_TN/TN_cov) > 0.1) or (Q30T_indel and (mut_reads_Q30T/Q30T_cov) > 0.1) or (TT_indel and (mut_reads_TT/TT_cov) > 0.1):
				continue
			elif ((mut_reads_TN-mut_reads_Q30N > 0.05*TN_cov) and (mut_reads_TN/TN_cov > 0.15)) or ((mut_reads_Q30N == 0) and (mut_reads_TN/TN_cov > 0.1)): #Filter positions with many mut reads in normal with low quality
				continue
			else: #We only accept substitutions
				norm_reads_number_TN, mut_reads_number_TN, normal_mapQ_TN, mut_mapQ_TN = prepare_mapQ(TN_bases_correct, TN_mapQ, alt_TT)
				norm_reads_number_TT, mut_reads_number_TT, normal_mapQ_TT, mut_mapQ_TT = prepare_mapQ(TT_bases_correct, TT_mapQ, alt_TT)
				normal_error_ratio, tumor_error_ratio = extract_error_ratio(line[4], mut_reads_number_TN, mut_reads_number_TT)
				normal_cigar, tumor_cigar = extract_cigar(line[5])
				max_dimer_value, GC_percentage = start(line[6], 'sequence')
				distance_position, mean_position = start(line[7], 'interval')
				mut_reads_Q30T = baseq_filter(samfile, line[8])
				if mut_reads_Q30T < args.tumor_mut_reads:
					continue
				else:
					print('%s_%s_%s' % (chrom, pos, args.name), Q30N_cov, Q30T_cov, mut_reads_Q30N, mut_reads_Q30T, normal_mapQ_TN, normal_mapQ_TT, mut_mapQ_TN, mut_mapQ_TT, normal_error_ratio, tumor_error_ratio, normal_cigar, tumor_cigar, max_dimer_value, GC_percentage, distance_position, mean_position, sep=',', end='\n')
		except ZeroDivisionError:
			continue

	Q30N.close()
	Q30T.close()
	TN.close()
	TT.close()
	stats.close()
	cigar.close()
	sequence.close()
	names.close()
