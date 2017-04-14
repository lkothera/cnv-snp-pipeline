#!/usr/bin/env python

# -----
# Title: pre_variant.py
# Author: John Phan (John.Phan@csra.com)
# -----

# Pre-Variant analysis before applying GATK as part of OTG-SNPCaller pipeline


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
		description='pre_variant.py'
	)
	parser.add_argument(
		'--in_dir',
		type=str,
		required=True,
		help='input directory'
	)
	parser.add_argument(
		'--out_dir',
		type=str,
		required=True,
		help='output directory'
	)
	parser.add_argument(
		'--base_dir',
		type=str,
		required=True,
		help='base directory'
	)
	args = parser.parse_args()
	return (
		args.in_dir,
		args.out_dir,
		args.base_dir
	)


def invoke(command):
	return Popen(command, stdout=PIPE, shell=True).stdout.read()


def run_pre_variant(
	in_dir, 
	out_dir,
	base_dir
):

	# get list of fastq files
	dirlist = glob.glob(in_dir+'/*rmDup_sorted.bam')

	mapping_stats = {}

	for bam_file in dirlist:
		basename = os.path.basename(bam_file)
		samplename = os.path.splitext(basename)[0]

		# run pre-variant script
		cmd = 'perl '+base_dir+'/scripts/proton_pre_snp_v1_AOS.pl -b '+bam_file+' -r 0.2 -pre '+out_dir+'/'+samplename
		sys.stderr.write(cmd+'\n')
		invoke(cmd)


def main():
	(
		in_dir,
		out_dir,
		base_dir
	) = parse_args()

	run_pre_variant(in_dir, out_dir, base_dir)

	return True


if __name__ == "__main__":
	main()


