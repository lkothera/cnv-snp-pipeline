#!/usr/bin/env python

# -----
# Title: bwa_align.py
# Author: John Phan (John.Phan@csra.com)
# -----

# bwa align and sort all fastq files in a directory


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
		description='bwa_align.py'
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
		'--ext',
		type=str,
		required=True,
		help='input extension'
	)
	parser.add_argument(
		'--database',
		type=str,
		required=True,
		help='bwa database path'
	)
	args = parser.parse_args()
	return (
		args.in_dir,
		args.out_dir,
		args.ext,
		args.database
	)


def invoke(command):
	return Popen(command, stdout=PIPE, shell=True).stdout.read()


def run_alignments(
	in_dir, 
	out_dir,
	ext,
	database
):

	# get list of fastq files
	dirlist = glob.glob(in_dir+'/*.'+ext)

	mapping_stats = {}

	for fastq_file in dirlist:
		basename = os.path.basename(fastq_file)
		samplename = os.path.splitext(basename)[0]

		# bwa align
		cmd = 'bwa mem -t 4 '+database+' '+fastq_file+' 2> /dev/null | samtools view -b - | samtools sort -O bam -T tmp > '+out_dir+'/'+samplename+'.bam'
		sys.stderr.write(cmd+'\n')
		invoke(cmd)

		# samstats
		cmd = 'samtools stats '+out_dir+'/'+samplename+'.bam | grep ^SN | cut -f 2- > '+out_dir+'/'+samplename+'.stats'
		invoke(cmd)

		# extract fields from samstats
		mapping_stats[samplename] = {}

		# -- sequences
		cmd = 'cat '+out_dir+'/'+samplename+'.stats | grep "^sequences:"'
		sequences = invoke(cmd)
		m = re.match(
			'.*:\s*(\d*)',
			sequences
		)
		if m:
			mapping_stats[samplename]['sequences'] = int(m.group(1))

		# -- reads mapped
		cmd = 'cat '+out_dir+'/'+samplename+'.stats | grep "^reads mapped:"'
		reads_mapped = invoke(cmd)
		m = re.match(
			'.*:\s*(\d*)',
			reads_mapped
		)
		if m:
			mapping_stats[samplename]['reads_mapped'] = int(m.group(1))

		# -- percent mapped
		mapping_stats[samplename]['perc_mapped'] = '%0.2f' % (100*mapping_stats[samplename]['reads_mapped']/float(mapping_stats[samplename]['sequences']))

		# -- avarage quality
		cmd = 'cat '+out_dir+'/'+samplename+'.stats | grep "^average quality:"'
		reads_mapped = invoke(cmd)
		m = re.match(
			'.*:\s*(.*)',
			reads_mapped
		)
		if m:
			mapping_stats[samplename]['average_quality'] = float(m.group(1))
		
		# -- average length
		cmd = 'cat '+out_dir+'/'+samplename+'.stats | grep "^average length:"'
		reads_mapped = invoke(cmd)
		m = re.match(
			'.*:\s*(\d*)',
			reads_mapped
		)
		if m:
			mapping_stats[samplename]['average_length'] = int(m.group(1))

	# print mapping stats
	print 'Sample\t# Sequences\t# Mapped Sequences\t% Mapped\tAverage Length\tAverage Quality'
	for sample in mapping_stats.keys():
		print (
			sample+'\t'+
			str(mapping_stats[sample]['sequences'])+'\t'+
			str(mapping_stats[sample]['reads_mapped'])+'\t'+
			str(mapping_stats[sample]['perc_mapped'])+'\t'+
			str(mapping_stats[sample]['average_length'])+'\t'+
			str(mapping_stats[sample]['average_quality'])
		)

def main():
	(
		in_dir,
		out_dir,
		ext,
		database
	) = parse_args()

	run_alignments(in_dir, out_dir, ext, database)

	return True


if __name__ == "__main__":
	main()


