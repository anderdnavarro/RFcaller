#!/usr/bin/env python3

import argparse
from sys import stdin, argv
from re import findall, match
from statistics import mean
from pysam import faidx
from regex import findall as regexfindall
from math import ceil

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Filter mutations extracted by samtools mpileup using a tumor and normal bam. The output are four files: two VCFs with the positions that have mutations (independent files for tumor and normal samples) and the tumor and normal pileup for these mutations. Be careful --> The pileup must not contain any extra field, only default ones!!')
	parser.add_argument(
	'-i', '--input', type=str, required=False, default= '-', metavar='MQ',
	help='Pileup file for tumor and normal in these order.')
	parser.add_argument(
	'-o', '--output_name', type=str, required=True, metavar='OUTPUT ID',
	help='Sample ID to appear in VCF')
	parser.add_argument(
	'-O', '--output_dir', type=str, required=True, metavar='OUTPUT DIR',
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
	'-Tmut', '--tumor_mut_reads', type=int, required=False, default=4, metavar='TUMOR MUT READS',
	help='Minimum number of mutated reads in tumor to select the position as mutated (default: %(default)s).')
	parser.add_argument(
	'-Nmut', '--normal_mut_reads', type=int, required=False, default=2, metavar='NORMAL MUT READS',
	help='Minimum number of mutated reads in normal to select the position as mutated (default: %(default)s).')
	parser.add_argument(
	'-c', '--contamination', type=float, required=False, default=0.05, metavar='CONTAMINATION',
	help='Percentage of tumor contamination in normal sample (default: %(default)s).')
	parser.add_argument(
	'-f', '--ref_path', type=str, required=True, metavar='REFERENCE GENOME',
	help='Path to the reference genome')
	parser.add_argument(
	'-r', '--region', type=int, default= 10, metavar='REGION',
	help='Window size to consider the region (default: %(default)s)')
	return parser.parse_args()

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
	return(set(indel_list), indel_list, variants)

def most_common_indel(single_indel_list, full_indel_list, cutoff, sample):
	alt, alt_count = [], []
	for indel in single_indel_list:
		count = full_indel_list.count(indel)
		if sample == 'tumor':
			if count >= cutoff:
				alt.append(indel)
				alt_count.append(str(count))
			else:
				continue
		elif sample == 'normal':
			if count > cutoff:
				alt.append(indel)
				alt_count.append(str(count))
			else:
				continue
	return(alt, alt_count)

def get_mut_position(variants, positions):
	positions_list = []
	for base in zip(variants, positions):
		if base[0] == ';': #To extract only those read names that contain the mutation we are looking for
			positions_list.append(int(base[1]))
		else:
			continue
	distance = max(positions_list)-min(positions_list)
	mean_distance = mean(positions_list)/args.reads_length
	return(distance, mean_distance)

def make_ref_list(alt_list, ref):
	ref_list = []
	for i in range(0, len(alt_list)):
		if len(alt_list[i]) == 1:
			ref_list.append(ref)
			continue
		else:
			nucleotides = ''.join(findall(r'[A-Z]', alt_list[i]))
			indel_deletion = match('-', alt_list[i])
			if indel_deletion is not None: #Deletions
				alt_list[i] = ref
				ref_list.append(ref + nucleotides)
				continue
			else: #Insertions
				alt_list[i] = ref + nucleotides
				ref_list.append(ref)
				continue
	return(ref_list, alt_list)

def check_indel(indel_bases):
	size_dict = {}

	## First we make a dictionary taking into account the length of the indel
	for indel in indel_bases:
		size_dict[len(indel)] = indel

	## Then we check if the smallest indel is inside the rest of the indels of the same position to consider only the smallest one
	n = 1
	smallest_indel = size_dict.get(min(size_dict))
	for key in size_dict:
		indel = size_dict.get(key)
		if smallest_indel not in indel:
			n += 1
			continue
		else:
			continue
	return (smallest_indel, n)

def decode_ASCII(code):
	decode_value = []
	for i in code:
		decode_value.append(ord(i)-33) # ASCII code -33 (samtools)
	return(float(mean(decode_value)))

def filter_indels(candidate_indels):
	for candidate in candidate_indels:
		column = vcf_dict.get(candidate)
		chrom, pos, output_name, ref_tum, alt_tum, cov_norm, BQN, MQN, indel_list_norm, cov_tum, BQT, MQT, alt_count_tum = str(column[0]), int(column[1]), str(column[2]), column[3].split(','), column[4].split(','), column[5], column[6], column[7], column[8].split(','), column[9], column[10], column[11], column[12].split(',')
		## Check wheter the mutation is only an indel and always the same type
		indel, deletion, insertion, first = True, False, False, True
		for base in zip(ref_tum, alt_tum):
			if len(base[0]) == 1 and len(base[1]) == 1:
				indel=False
				break
			elif len(base[0]) > 1 and len(base[1]) == 1 and first:
				deletion, first = True, False
				continue
			elif len(base[0]) == 1 and len(base[1]) > 1 and first:
				insertion, first = True, False
				continue
			elif len(base[0]) > 1 and len(base[1]) == 1 and insertion:
				indel=False
				break
			elif len(base[0]) == 1 and len(base[1]) > 1 and deletion:
				indel=False
				break
			else:
				continue
		
		## If it is always an insertion or a deletion and there isn't any indel in the same position in the normal
		if indel and ND_check_vcf_dict.get(candidate) == None:
			if deletion:
				ref, n  = check_indel(ref_tum)
				index_tum = ref_tum.index(ref)
				alt = alt_tum[0]
				real_indel = ref
				alt_count_norm = indel_list_norm.count('-' + str(len(real_indel)-1) + real_indel[1:])
			elif insertion:
				ref = ref_tum[0]
				alt, n = check_indel(alt_tum)
				index_tum = alt_tum.index(alt)
				real_indel = alt
				alt_count_norm = indel_list_norm.count('+' + str(len(real_indel)-1) + real_indel[1:])
			else:
				continue

			## Calculate MAF
			alt_count_tum = int(alt_count_tum[index_tum])
			MAF = alt_count_tum/cov_tum

			## If there are different indels for the same position in normal, it is removed
			try:
				m = sum(1 for x in range(pos-args.region, pos+args.region) if x in ND_check_dict.get(chrom))
			except TypeError:
				m = 0
				
			if n != 1 or m != 0:
				continue
			else:
				# Write mini.pileup
				ND_pileup.write('\t'.join(ND_pileup_dict.get(candidate))+'\n')
				TD_pileup.write('\t'.join(TD_pileup_dict.get(candidate))+'\n')
				# Write interval size
				mutations_interval.write('\t'.join(interval_dict.get(candidate))+'\n')
				# Print mini.pileup.vcf
				print(chrom, pos, output_name, ref, alt, '.', '.', 'VAF=%.3f' % float(MAF), 'DP:BQ:MQ:AD', '%i:%.1f:%.1f:%i' % (cov_norm, BQN, MQN, alt_count_norm), '%i:%.1f:%.1f:%i' % (cov_tum, BQT, MQT, alt_count_tum), sep='\t', end='\n')

				## Create the indel_features file
				# Is the indel repeat in the region?
				if len(real_indel) == 2:
					repeat_indel = 0
				else:
					region = chrom + ':' + str(pos-15) + '-' + str(pos+15)
					sequence = faidx(args.ref_path, region).splitlines() #To search the sequence in the reference genome
					repeat_indel = len(regexfindall(real_indel[1:], sequence[1], overlapped=True))

				# Write the features to output
				features.write('%s\t%s\n' % (candidate+'_'+args.output_name, repeat_indel))
		else:
			continue

def main_function(pileup):
	for line in pileup:
		column = line.strip().split()
		chrom, pos, ref, cov_tum, variants_tum, baseq_tum, positions_tum, mapping_tum, cov_norm, variants_norm, baseq_norm, positions_norm, mapping_norm = str(column[0]), column[1], column[2], int(column[3]), column[4].upper(), column[5], column[6].split(','), column[7], int(column[8]), column[9].upper(), column[10], column[11].split(','), column[12]
		key_id = chrom + '_' + str(pos)

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
			single_indel_tum, indel_list_tum, new_variants_tum = correct_variants_list(variants_tum)
			single_indel_norm, indel_list_norm, new_variants_norm = correct_variants_list(variants_norm)
			## Filter positions with no real mutations in tumor
			if len(single_indel_tum) == 0:
				continue
			else:
				alt_tum, alt_count_tum = most_common_indel(single_indel_tum, indel_list_tum, args.tumor_mut_reads, 'tumor')
				if len(alt_count_tum) == 0:
					continue
				else:
					# Calculate tumor's MAF
					MAF_tum = int(alt_count_tum[0])/cov_tum
					pass
			
			## Check the number of variants for normal
			if len(single_indel_norm) == 0: # If there isn't any mutation in the normal
				alt_norm, alt_count_norm = [], ['0']
				pass
			else:
				## Updates the number of mut reads allowed in the normal taking into account the contamination
				max_normal_mut_reads = ceil(args.normal_mut_reads + (cov_norm*args.contamination*MAF_tum))
				# Filter
				alt_norm, alt_count_norm = most_common_indel(single_indel_norm, indel_list_norm, max_normal_mut_reads, 'normal')
				if len(alt_norm) != 0: # If there is any mutation in the normal with enough reads, save the position
					ref_norm, alt_norm = make_ref_list(alt_norm, ref)
					ND_check_vcf_dict[key_id] = [ref_norm, alt_norm]
					try:
						ND_check_dict[chrom].append(int(pos))
					except KeyError:
						ND_check_dict[chrom] = [int(pos)]
				else:
					pass
			
			## Save the information to print
			if len(alt_tum) != 0:
				TD_list.append(key_id)
				## Calculate the BQ and MQ for the position of the mutation
				BQN, BQT = decode_ASCII(baseq_norm), decode_ASCII(baseq_tum)
				MQN, MQT = decode_ASCII(mapping_norm), decode_ASCII(mapping_tum)
				# Save the pileups
				TD_pileup_dict[key_id] = [str(x) for x in [chrom, pos, ref, cov_tum, variants_tum, baseq_tum, mapping_tum, ','.join(positions_tum)]]
				ND_pileup_dict[key_id] = [str(x) for x in [chrom, pos, ref, cov_norm, variants_norm, baseq_norm, mapping_norm, ','.join(positions_norm)]]
				# Save the VCF
				ref_tum, alt_tum = make_ref_list(alt_tum, ref)
				vcf_dict[key_id] = [chrom, pos, args.output_name, ','.join(ref_tum), ','.join(alt_tum), cov_norm, BQN, MQN, ','.join(indel_list_norm), cov_tum, BQT, MQT, ','.join(alt_count_tum)]
				# Check in which position of the reads the mutation appears
				distance_position, mean_position = get_mut_position(new_variants_tum, positions_tum)
				interval_dict[key_id] = [str(x) for x in [key_id+'_'+args.output_name, distance_position, mean_position]]
			else:
				continue
		else:
			continue

if __name__ == '__main__':
	## Arguments
	args = parse_args()

	## Create empty variables
	vcf_dict, TD_pileup_dict, ND_pileup_dict, interval_dict, ND_check_dict, ND_check_vcf_dict, ND_count_dict = {}, {}, {}, {}, {}, {}, {}
	TD_list = []

	## Outputs
	ND_pileup = open(args.output_dir+'/'+args.ND_name+'.mini.pileup', 'w+')
	TD_pileup = open(args.output_dir+'/'+args.TD_name+'.mini.pileup', 'w+')
	mutations_interval = open(args.output_dir+'/'+args.output_name+'.mutations.interval', 'w+')
	features = open(args.output_dir+'/'+args.output_name+'.features', 'w+')
	

	## Headers
	# VCF
	arg_list = []
	for arg in sorted(vars(args)):
		arg_list.append('--'+arg)
		arg_list.append(str(getattr(args, arg)))
	header = ['#CHROM','POS','ID','REF','ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.ND_name, args.TD_name]
	print('##fileformat=VCFv4.2', '##INFO=<ID=VAF,Number=1,Type=Float,Description="Allele frequency of the alternative allele">', '##FILTER=<ID=DPT_INDEL,Description="Depth coverage in tumor >=%i for INDELs">' % args.tumor_coverage, '##FILTER=<ID=DPN_INDEL,Description="Depth coverage in normal >=%i for INDELs">' % args.normal_coverage, '##FILTER=<ID=ADT_INDEL,Description="Minimum number of alternative reads in tumor >=%i for INDELs">' % args.tumor_mut_reads, '##FILTER=<ID=ADN_INDEL,Description="Maximum number of alternative reads in normal <=%i for INDELs">' % args.normal_mut_reads, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth coverage">',	'##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average base quality">', '##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Average mapping quality">', '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Number of alternative reads">', '##Command=python3 %s %s' % (argv[0], ' '.join(arg_list)), '\t'.join(header), sep='\n')
	# Interval
	mutations_interval.write('ID' + "\t" + 'INTERVAL_SIZE' + "\t" + 'MEAN_POSITION' + "\n")
	# Features
	features.write('ID' + "\t" + 'N_repeat_indel' + "\n")

	## Take the input
	if args.input == '-':
		with stdin as pileup:
			main_function(pileup)
	else:
		with open(args.input) as pileup:
			main_function(pileup)

	## Filter the mini.pileup
	filter_indels(TD_list)

	## Close output files
	ND_pileup.close()
	TD_pileup.close()
	mutations_interval.close()
	features.close()

