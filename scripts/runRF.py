#!/usr/bin/env python3

import argparse
import pandas as pd
from joblib import load

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script launches a Machine Learning algorithm to detect Somatic Mutations')
	parser.add_argument(
	'-i', '--input', type=str, required=True, metavar='CSV',
	help='Dataset in CSV format')
	parser.add_argument(
	'-a', '--algorithm', type=str, required=True, metavar='ALGORITHM',
	help='The trained algorithm to detect somatic mutations')
	return parser.parse_args()

if __name__ == '__main__':

	## Arguments
	args = parse_args()

	## Load the trained algorithm
	model = load(args.algorithm)

	## Input CSV data
	dataset = pd.read_csv(args.input, sep=',')
	# Extract the ID of each mutation to understand the results
	names = list(dataset['#ID'].values)
	del dataset['#ID']

	## Define the dataset
	X = dataset.values
	
	## Use the algorithm
	algorithm_result = list(model.predict(X))

	## Print the output
	for case in zip(names, algorithm_result):
		print(case[0], case[1], sep='\t', end='\n')
