#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from joblib import dump
from sklearn.metrics import mean_absolute_error
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import MaxAbsScaler, StandardScaler
from tpot.builtins import StackingEstimator
from tpot.export_utils import set_param_recursive

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='This script trains a regression machine learning algorithm to detect INDELs')
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

# Average CV score on the training set was: -0.07972758187757621
exported_pipeline = make_pipeline(
	StackingEstimator(estimator=GradientBoostingRegressor(alpha=0.8, learning_rate=0.001, loss="lad", max_depth=3, max_features=0.6000000000000001, min_samples_leaf=16, min_samples_split=19, n_estimators=100, subsample=0.3)),
	StackingEstimator(estimator=RidgeCV()),
	StackingEstimator(estimator=KNeighborsRegressor(n_neighbors=28, p=2, weights="uniform")),
	StandardScaler(),
	MaxAbsScaler(),
	RandomForestRegressor(bootstrap=False, max_features=0.25, min_samples_leaf=3, min_samples_split=3, n_estimators=100)
)

## Fix random state for all the steps in exported pipeline
exported_pipeline.fit(training_features, training_target)

## Export algorithm
dump(exported_pipeline, args.output)

## Scores
mse_test = mean_absolute_error(testing_target, exported_pipeline.predict(testing_features))
print("Test Set Mean Absolute Error: %.5f" % mse_test)
