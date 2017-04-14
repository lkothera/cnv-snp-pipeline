#!/usr/bin/env python

# -----
# Title: indel_realign.py
# Author: John Phan (John.Phan@csra.com)
# -----

# GATK Indel Realignment


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
		description='tmap_align.py'
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
		'--reference',
		type=str,
		required=True,
		help='reference sequence path'
	)
	args = parser.parse_args()
	return (
		args.in_dir,
		args.out_dir,
		args.reference
	)


def invoke(command):
	return Popen(command, stdout=PIPE, shell=True).stdout.read()


def run_indel_realignments(
	in_dir, 
	out_dir,
	reference
):

	# get list of fastq files
	dirlist = glob.glob(in_dir+'/*rmDup_sorted.bam')

	mapping_stats = {}

	for bam_file in dirlist:
		basename = os.path.basename(bam_file)
		samplename = os.path.splitext(basename)[0]

		# GATK Realigner Target Creator
		cmd = 'GenomeAnalysisTK.sh -T RealignerTargetCreator -l INFO -I '+bam_file+' -R '+reference+' -o '+out_dir+'/'+samplename+'.intervals'
		sys.stderr.write(cmd+'\n')
		invoke(cmd)

		# GATK Indel Realigner
		cmd = 'GenomeAnalysisTK.sh -T IndelRealigner -l INFO -I '+bam_file+' -R '+reference+' -targetIntervals '+out_dir+'/'+samplename+'.intervals -o '+out_dir+'/'+samplename+'.sort.fix.bam'
		sys.stderr.write(cmd+'\n')
		invoke(cmd)

		# Update read groups
		cmd = 'picard AddOrReplaceReadGroups I='+out_dir+'/'+samplename+'.sort.fix.bam O='+out_dir+'/'+samplename+'.sort.fix.rg.bam LB=lib PL=illumina PU=barcode SM='+samplename
		sys.stderr.write(cmd+'\n')
		invoke(cmd)

		# Index
		cmd = 'samtools index '+out_dir+'/'+samplename+'.sort.fix.rg.bam'
		sys.stderr.write(cmd+'\n')
		invoke(cmd)



def main():
	(
		in_dir,
		out_dir,
		reference
	) = parse_args()

	run_indel_realignments(in_dir, out_dir, reference)

	return True


if __name__ == "__main__":
	main()


