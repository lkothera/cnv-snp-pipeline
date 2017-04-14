#!/usr/bin/env python

# -----
# Title: combine_exports.py
# Author: John Phan (John.Phan@csra.com)
# -----

# combine exported bed files from conifer


import os
import glob
import re
import sys
import argparse
from pprint import pprint
from subprocess import (PIPE, Popen)


def parse_args():
	# parse command line arguments and return
	# as list of argument values
	parser = argparse.ArgumentParser(
		description='combine_exports.py'
	)
	parser.add_argument(
		'--in_dir',
		type=str,
		required=True,
		help='input directory'
	)
	parser.add_argument(
		'--out',
		type=str,
		required=True,
		help='output data file'
	)
	args = parser.parse_args()
	return (
		args.in_dir,
		args.out
	)


def invoke(command):
	return Popen(command, stdout=PIPE, shell=True).stdout.read()


def read_bed_files(
	in_dir 
):

	bed_data = {}
	probe_list = {}

	# get list of fastq files
	dirlist = glob.glob(in_dir+'/*.bed')

	for bed_file_name in dirlist:
		basename = os.path.basename(bed_file_name)
		sample_name = os.path.splitext(basename)[0]
		sample_name = os.path.splitext(sample_name)[0]

		# open bed file
		try:
			bed_file = open(bed_file_name, 'rU')
		except IOError as err:
			print (
				'!error: '
				'read_bed_files(): '
				'cannot open bed file: '+bed_file_name
			)
			print err
			return False


		try:
			while True:
				line = bed_file.next().strip()
				line_elems = re.split(r'\t', line)

				probe_start = line_elems[1]
				gene_name = line_elems[3]
				cnv_value = float(line_elems[4])
				probe_name = gene_name+'_'+probe_start
				probe_list[probe_name] = 1

				if sample_name not in bed_data:
					bed_data[sample_name] = {}

				bed_data[sample_name][probe_name] = cnv_value

		except(StopIteration):
			bed_file.close

	return (bed_data, probe_list)


def main():
	(
		in_dir,
		out
	) = parse_args()

	(data, probe_names) = read_bed_files(in_dir)


	try:
		out_file = open(out, 'w+')
	except IOError as err:
		print (
			'!error: '
			'combine_exports.py: '
			'cannot open output file: '+out
		)
		print err
		sys.exit(1)

	# print header
	for sample in sorted(data):
		out_file.write('\t'+sample)
	out_file.write('\n')

	for probe in sorted(probe_names):
		out_file.write(probe)
		for sample in sorted(data):
			out_file.write('\t'+str(data[sample][probe]))
		out_file.write('\n')

	out_file.close()

	sys.exit(0)


if __name__ == "__main__":
	main()


