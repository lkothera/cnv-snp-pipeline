#!/usr/bin/env python

# -----
# Title: extract_label_samples.py
# Author: John Phan (John.Phan@csra.com)
# -----

# extract and label samples from a tab-delimited text file


import os
import re
import sys
import argparse
from pprint import pprint


def parse_args():

	# parse command line arguments and return
	# as list of argument values
	parser = argparse.ArgumentParser(
		description='extract_label_samples.py'
	)
	parser.add_argument(
		'--data_path',
		type=str,
		required=True,
		help='input sample data file'
	)
	parser.add_argument(
		'--label_path',
		type=str,
		required=True,
		help='sample to label file'
	)
	parser.add_argument(
		'--column',
		type=int,
		required=True,
		help='column in sample to label file'
	)
	args = parser.parse_args()
	return (
		args.data_path,
		args.label_path,
		args.column
	)


def load_sample_to_label_map(label_path, column):

	sample_to_label_map = {}

	# open label file
	try:
		label_file = open(label_path, 'rU')
	except IOError as err:
		print (
			'!error: '
			'load_sample_to_label_map(): '
			'cannot open label file: '+label_path
		)
		print err
		return False

	# read first header line
	header = label_file.next()

	# iterate through remaining lines
	try:
		while True:
			line = label_file.next().strip()
			line_elems = re.split(r'\t', line)
			sample_name = line_elems[0]
			file_present = int(line_elems[4])
			label = int(line_elems[4+column])
			if file_present > 0 and label > 0:
				sample_to_label_map[sample_name] = label
	except(StopIteration):
		label_file.close()

	return sample_to_label_map


def load_data(data_path):

	data = {}
	sample_columns = {}

	# open data file
	try:
		data_file = open(data_path, 'rU')
	except IOError as err:
		print (
			'!error: '
			'load_data(): '
			'cannot open data file: '+data_path
		)
		print err
		return False
	
	# read the header line
	header = data_file.next().strip()
	
	# parse sample columns, keep track of the column to which each sample belongs
	header_elems = re.split(r'\t', header)
	for i in range(0, len(header_elems)):
		sample_columns[header_elems[i]] = i

	# iterate through remaining lines
	try:
		while True:
			line = data_file.next().strip()
			line_elems = re.split(r'\t', line)
			feature_name = line_elems[0]
			data[feature_name] = line_elems[1:]
	except(StopIteration):
		data_file.close()

	return (
		sample_columns,
		data
	)


def print_labeled_data(
	sample_columns,
	sample_to_label_map,
	data
):

	# create list of samples to include
	sample_list = []
	for sample in sorted(sample_to_label_map):
		if sample in sample_columns:
			sample_list.append(sample)
	
	# print header
	for sample in sample_list:
		sys.stdout.write('\t'+sample)
	sys.stdout.write('\n')

	# print labels
	sys.stdout.write('class')
	for sample in sample_list:
		sys.stdout.write('\t'+str(sample_to_label_map[sample]))
	sys.stdout.write('\n')

	# print each feature
	for feature in sorted(data):
		sys.stdout.write(feature)
		for sample in sample_list:
			sys.stdout.write('\t'+data[feature][sample_columns[sample]])
		sys.stdout.write('\n')

	return True


def main():

	(
		data_path,
		label_path,
		column
	) = parse_args()

	sample_to_label_map = load_sample_to_label_map(label_path, column)
	if not sample_to_label_map:
		print (
			'!error: '
			'extract_label_samples.py: '
			'cannot load sample to label map'
		)
		sys.exit(1)

	#pprint(sample_to_label_map)

	(sample_columns, data) = load_data(data_path)
	if not data:
		print (
			'!error: '
			'extract_label_samples.py: '
			'cannot load sample to label map'
		)
		sys.exit(1)

	#pprint(sample_columns)
	#pprint(data)

	if not print_labeled_data(
		sample_columns,
		sample_to_label_map,
		data
	):
		print (
			'!error: '
			'extract_label_samples.py: '
			'cannot print labeled data'
		)
		sys.exit(1)

	sys.exit(0)


if __name__ == "__main__":
	main()

