#!/usr/bin/env python

# -----
# Title: assoc_variants.py
# Author: John Phan (John.Phan@csra.com)
# Last Modified: 2016/10/03
# -----

# associate variants in vcf to groups


import os
import glob
import sys
import argparse
from pprint import pprint
from subprocess import (PIPE, Popen)
import vcf
import scipy.stats as stats

def parse_args():
	# parse command line arguments and return
	# as list of argument values
	parser = argparse.ArgumentParser(
		description='assoc_variants.py'
	)
	parser.add_argument(
		'--vcf',
		type=str,
		required=True,
		help='VCF file'
	)
	parser.add_argument(
		'--compare',
		type=str,
		required=True,
		help='path to file that contains groupings'
	)
	parser.add_argument(
		'--column',
		type=int,
		required=True,
		help='column in compare file to extract groupings'
	)
	args = parser.parse_args()
	return (
		args.vcf,
		args.compare,
		args.column
	)


def invoke(command):
	return Popen(command, stdout=PIPE, shell=True).stdout.read()


def read_compare(
	compare_path,
	column
):
	compare_dict = {}

	try:
		compare_file = open(compare_path, 'rU')
	except IOError as err:
		print (
			"!error: "
			"read_compare(): "
			"cannot read compare file: "+compare_path
		)
		print err
		return False

	try:
		line = compare_file.next()	# skip header line
		while True:
			line = compare_file.next().strip()

			if line == None:
				continue

			line_elems = line.split('\t')

			# column 1 corresponds to 0-indexed column 5
			col = column+4

			# get sample name and group
			sample = line_elems[0]
			group = line_elems[col]

			compare_dict[sample] = int(group)
	except(StopIteration):
		compare_file.close()

	return compare_dict	


def main():
	(
		vcf_path,
		compare,
		column
	) = parse_args()

	compare_dict = read_compare(compare, column)
	if not compare_dict:
		print (
			"!error: "
			"assoc_variants.py: "
			"read_compare() failed"
		)
		sys.exit(1)

	#pprint(compare_dict)

	print "ID\tPOS\tREF\tALT\t# Called\t# Heterozygous\t# Homozygous Alt. Allele\t# Homozygous Ref. Allele\t# Ref. Susceptible (Group 1)\t# Ref. Resistant (Group 2)\t# Alt. Susceptible (Group 1)\t# Alt. Resistant (Group 2)\tP-Val., Fisher Exact\tFE Statistic\tOdds Ratio\tAlt. Susceptible Samples\tAlt. Resistant Samples"

	vcf_reader = vcf.Reader(open(vcf_path, 'r'));
	for record in vcf_reader:

		#print record.REF
		#print record.CHROM
		res_ref = 0
		res_alt = 0
		sus_ref = 0
		sus_alt = 0
		num_het = 0
		num_hom_ref = 0
		num_hom_alt = 0
		sus_alt_list = []
		res_alt_list = []

		for call in record.samples:
			# extract sample name
			sample_name = call.sample.split(".")[0]
			#print(sample_name)
			if compare_dict[sample_name] == 0:
				# sample not included in comparison
				continue

			if call.gt_type is None:
				# no call was made on the sample
				continue 

			if call.gt_type == 0:
				num_hom_ref += 1
				# reference type
				if compare_dict[sample_name] == 1:
					# susceptible
					sus_ref += 1
				else:
					# resistant
					res_ref += 1
			else:
				if call.gt_type == 1:
					num_het += 1
				else:
					num_hom_alt += 1

				# alt type
				if compare_dict[sample_name] == 1:
					# susceptible
					sus_alt += 1
					sus_alt_list.append(sample_name)
				else:
					res_alt += 1
					res_alt_list.append(sample_name)

		#print "Contingency: [" + str(res_ref) + ", " + str(sus_ref) + "], [" + str(res_alt) + ", " + str(sus_alt) + "]"
		oddsratio, pvalue = stats.fisher_exact([[res_ref, sus_ref], [res_alt, sus_alt]])
		#print pvalue

		stat = 0
		try:
			stat = res_ref/float(res_ref+sus_ref)-res_alt/float(res_alt+sus_alt)
		except ZeroDivisionError:
			# p-value will be non-significant anyway
			stat = 0

		print (
			record.CHROM + "\t" +
			str(record.POS) + "\t" +
			record.REF + "\t" +
			",".join(map(str, record.ALT)) + "\t" +
			str(num_het+num_hom_ref+num_hom_alt) + "\t" +
			str(num_het) + "\t" +
			str(num_hom_alt) + "\t" +
			str(num_hom_ref) + "\t" +
			str(sus_ref) + "\t" +
			str(res_ref) + "\t" +
			str(sus_alt) + "\t" +
			str(res_alt) + "\t" +
			str(pvalue) + "\t" +
			str(stat) + "\t" +
			str(oddsratio) + "\t" +
			",".join(sus_alt_list) + "\t" +
			",".join(res_alt_list)
		)


if __name__ == "__main__":
	main()


