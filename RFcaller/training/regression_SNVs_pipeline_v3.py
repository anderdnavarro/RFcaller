#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from copy import copy
from joblib import dump
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.metrics import mean_absolute_error
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor
from sklearn.svm import LinearSVR
from tpot.builtins import StackingEstimator
from sklearn.preprocessing import FunctionTransformer

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script trains a regression machine learning algorithm to detect SNVs')
	parser.add_argument(
	'-i', '--input', type=str, required=True, metavar='CSV',
	help='A training dataset in CSV format. Make sure that the outcome column is labeled "Result" in the data file')
	parser.add_argument(
	'-o', '--output', type=str, required=True, metavar='OUTPUT',
	help='Output algorithm')
	return parser.parse_args()

## Arguments
args = parse_args()

## Input file
tpot_data = pd.read_csv(args.input, sep=',')
# Extract some columns
tpot_data_results = tpot_data['Result'].values
names = tpot_data['#ID'].values
# Remove unnecessary data
del tpot_data['Result']
del tpot_data['#ID']
features = tpot_data.values

## Training set
training_features, testing_features, training_target, testing_target = train_test_split(features, tpot_data_results, random_state=40)

## Average CV score on the training set was:-0.0061893517870134935
exported_pipeline = make_pipeline(
	make_union(
		make_union(
			FunctionTransformer(copy),
			make_union(
				make_union(
					make_union(
						StackingEstimator(estimator=GradientBoostingRegressor(alpha=0.95, learning_rate=0.5, loss="ls", max_depth=1, max_features=0.45, min_samples_leaf=6, min_samples_split=14, n_estimators=100, subsample=0.6500000000000001)),
						FunctionTransformer(copy)
					),
					StackingEstimator(estimator=LinearSVR(C=0.5, dual=False, epsilon=0.1, loss="squared_epsilon_insensitive", tol=1e-05))
				),
				FunctionTransformer(copy)
			)
		),
		FunctionTransformer(copy)
	),
	ExtraTreesRegressor(bootstrap=False, max_features=0.6500000000000001, min_samples_leaf=2, min_samples_split=9, n_estimators=100)
)

## Fix random state for all the steps in exported pipeline
exported_pipeline.fit(training_features, training_target)

## Export algorithm
dump(exported_pipeline, args.output)

## Scores
mse_test = mean_absolute_error(testing_target, exported_pipeline.predict(testing_features))
print("Test Set Mean Absolute Error: %.5f" % mse_test)
