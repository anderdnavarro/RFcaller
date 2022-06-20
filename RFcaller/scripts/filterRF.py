#!/usr/bin/env python3

from sys import argv
import argparse
from math import sqrt, ceil, floor

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script parses the result from the machine learning algorithm')
	parser.add_argument(
	'-i', '--input', type=str, required=True, metavar='INPUT',
	help='Input file containing the result from the machine learning algorithm')
	parser.add_argument(
	'-p', '--pileup', type=str, required=True, metavar='PILEUP_VCF',
	help='The mini.pileup VCF')
	parser.add_argument(
	'-v', '--vcf', type=str, required=True, metavar='VCF',
	help='The bcftools VCF')
	parser.add_argument(
	'-t', '--threshold', type=float, required=True, metavar='THRESHOLD',
	help='Minimum regression value to consider a mutation as good (default: %(default)s)')
	parser.add_argument(
	'--pipeline', type=str, required=True, choices=['SNV', 'INDEL'], metavar='PIPELINE',
	help='Which pipeline is running')
	parser.add_argument(
	'--polyIndel_threshold', type=float, required=False, default = 0, metavar='POLYINDEL_THRESHOLD',
	help='Regression value to consider a polyIndel as good')
	parser.add_argument(
	'--polyIndel_list', type=str, required=False, default = None, metavar='POLYINDEL_LIST',
	help='BED list with polyIndel positions')
	parser.add_argument(
	'-c', '--contamination', type=float, required=False, default=0.05, metavar='CONTAMINATION',
	help='Percentage of tumor contamination in normal sample (default: %(default)s).')
	return parser.parse_args()

def check_polyIndel(polyList, key):
	if len(polyList) == 0:
		return(False)
	else:
		if key in polyList:
			return(True)
		else:
			return(False)

def main_function(result, pileup, vcf):
	for line in zip(result, pileup, vcf):
		pileup_col, vcf_col =  line[1].strip(), line[2].strip()
		## Process the header of the pileup VCF and add new lines 
		while pileup_col.startswith('#'):
			if pileup_col.startswith('##'):
				if pileup_col.startswith('##fileformat'):
					print(pileup_col, '##FILTER=<ID=PASS,Description="All filters passed">', '##FILTER=<ID=LIKELY_GERMINAL,Description="Excess contamination in the normal">', '##INFO=<ID=RF,Number=1,Type=Float,Description="Regression value from the random forest algorithm">', sep='\n')
				else:
					print(pileup_col)
			else:
				print('##Command=python3 %s %s' % (argv[0], ' '.join(arg_list)), pileup_col, sep='\n')
			pileup_col = next(pileup).strip()

		## Process the header of the bcftools VCF
		while vcf_col.startswith('#'):
			vcf_col = next(vcf).strip()
		
		## Split lines
		res_col, pileup_col, vcf_col = line[0].strip().split(), pileup_col.split(), vcf_col.split()
		res_id, res_value, pileup_id, vcf_id, vcf_qual = '_'.join(res_col[0].split('_')[:2]), float(res_col[1]), '_'.join(pileup_col[:2]), '_'.join(vcf_col[:2]), float(vcf_col[5])
		
		## Check that there is the same mutation in both files
		while res_id != pileup_id or res_id != vcf_id:
			if res_id != pileup_id:
				pileup_col = next(pileup).strip().split()
				pileup_id = '_'.join(pileup_col[:2])
			elif res_id != vcf_id:
				vcf_col = next(vcf).strip().split()
				vcf_id, vcf_qual = '_'.join(vcf_col[:2]), float(vcf_col[5])
				
		## Cutoff
		if args.pipeline == 'SNV':
			real_qual = vcf_qual*(res_value**2)
			threshold = args.threshold
		elif args.pipeline == 'INDEL':
			polyIndel = check_polyIndel(polyList, res_id)
			if polyIndel == True:
				real_qual = res_value
				threshold = args.polyIndel_threshold
			else:
				real_qual = vcf_qual**res_value
				threshold = args.threshold
		
		if real_qual >= threshold:
			tumor_maf, normal_cov, normal_mut, tumor_mut = float(pileup_col[7].split('=')[1]), int(pileup_col[9].split(':')[0]), int(pileup_col[9].split(':')[3]), int(pileup_col[10].split(':')[3])
			if args.contamination != 0:
				expected_normal_maf = float(tumor_maf*args.contamination)
				cte = 1.96*sqrt(expected_normal_maf*(1-expected_normal_maf)/normal_cov)
				if normal_cov > 10:
					max_mut_reads_normal = int(ceil((expected_normal_maf+cte)*normal_cov))
				else:
					max_mut_reads_normal = int(floor((expected_normal_maf+cte)*normal_cov))
			else:
				max_mut_reads_normal = 1
				
			if (normal_mut > max_mut_reads_normal) or (tumor_mut <= 5 and normal_mut >= 2):
				print('\t'.join(pileup_col[:5]), '%.4f' % real_qual, 'LIKELY_GERMINAL', '%s;RF=%.4f' % (pileup_col[7], res_value), '\t'.join(pileup_col[8:]), sep='\t')
			else:
				print('\t'.join(pileup_col[:5]), '%.4f' % real_qual, 'PASS', '%s;RF=%.4f' % (pileup_col[7], res_value), '\t'.join(pileup_col[8:]), sep='\t')
		else:
			continue

if __name__ == '__main__':

	## Arguments
	args = parse_args()

	## New VCF command
	arg_list = []
	for arg in sorted(vars(args)):
		arg_list.append('--'+arg)
		arg_list.append(str(getattr(args, arg)))

	## Save polyIndel list if exists
	polyList = []
	if args.polyIndel_list is None:
		pass
	else:
		with open(args.polyIndel_list) as polyFile:
			for line in polyFile:
				line = line.strip().split()
				polyList.append(str(line[0])+'_'+str(line[1]))

	## Open the files
	with open(args.input) as result, open(args.pileup) as pileup, open(args.vcf) as vcf:
		main_function(result, pileup, vcf)
