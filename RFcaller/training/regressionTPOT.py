#!/usr/bin/env python3

import argparse
import pandas as pd
from tpot import TPOTRegressor
from sklearn.model_selection import train_test_split

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script selects the best pipeline to train a regression algorithm')
	parser.add_argument(
	'-i', '--input', type=str, required=True, metavar='CSV',
	help='A training dataset in CSV format. Make sure that the outcome column is labeled "Result" in the data file')
	parser.add_argument(
	'-o', '--output', type=str, required=True, metavar='OUTPUT',
	help='Output pipeline')
	parser.add_argument(
	'-t', '--threads', type=int, required=False, default=10, metavar='THREADS',
	help='Number of threads to use (default: %(default)s).')
	return parser.parse_args()

## Arguments
args = parse_args()

## Input data
dataset = pd.read_csv(args.input, sep=',')
# Get the solution
Y = dataset['Result'].values
# Remove unnecessary data
del dataset['Result']
del dataset['#ID']
X = dataset.values

## Random split of data into train_data and test_data
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, shuffle=True)

## Create the function which is going to optimize the regression pipeline
pipeline_optimizer = TPOTRegressor(generations=150, population_size=150, cv=5, random_state=42, n_jobs=args.threads, verbosity=2)

## Search for the best regrssion pipeline
pipeline_optimizer.fit(X_train, Y_train)

## Print the score of the pipeline
print(pipeline_optimizer.score(X_test, Y_test))

## Export the pipeline
pipeline_optimizer.export(args.output)
