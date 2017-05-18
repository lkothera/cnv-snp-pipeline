#!/usr/bin/env python

# -----
# Title: generate_single_contig.py
# Author: John Phan (John.Phan@csra.com)
# -----

# combine all exons into a single contig and generate a probe file for conifer


import os
import glob
import re
import sys
import argparse
from pprint import pprint
from subprocess import (PIPE, Popen)
from Bio import SeqIO

def parse_args():
	# parse command line arguments and return
	# as list of argument values
	parser = argparse.ArgumentParser(
		description='generate_single_contig.py'
	)
	parser.add_argument(
		'--in_fasta',
		type=str,
		required=True,
		help='input fasta file'
	)
	parser.add_argument(
		'--out_fasta',
		type=str,
		required=True,
		help='output fasta file'
	)
	parser.add_argument(
		'--in_probe',
		type=str,
		required=True,
		help='input probe file'
	)
	parser.add_argument(
		'--out_probe',
		type=str,
		required=True,
		help='output probe file'
	)
	parser.add_argument(
		'--out_gene',
		type=str,
		required=True,
		help='output gene coord file'
	)
	args = parser.parse_args()
	return (
		args.in_fasta,
		args.out_fasta,
		args.in_probe,
		args.out_probe,
		args.out_gene
	)


def invoke(command):
	return Popen(command, stdout=PIPE, shell=True).stdout.read()


def load_probes(
	in_probe
):

	# open input probe file and read into a data structure
	try:
		in_probe_file = open(in_probe, 'rU')
	except IOError as err:
		print (
			'!error: '
			'load_probes(): '
			'cannot open input probe file: '+in_probe
		)
		print err
		return False
	
	# first line is header
	line = in_probe_file.next()

	# initialize dict data structure
	super_contigs = {}
		
	try:

		while True:
			line = in_probe_file.next().strip()
			line_elems = re.split(r'\t', line)
			
			contig_name = line_elems[0]
			m = re.match(
				'supercont.\.(.*)',
				contig_name
			)
			contig_number = 0
			if m:
				contig_number = int(m.group(1))
			else:
				print (
					'!warning: '
					'load_probes(): '
					'invalid supercontig name in probe file:'+line_elems[0]
				)
				continue
			start = int(line_elems[1])
			stop = int(line_elems[2])
			gene = line_elems[3]

			if contig_number not in super_contigs:
				super_contigs[contig_number] = {}
				super_contigs[contig_number]['orig_name'] = contig_name
				super_contigs[contig_number]['probes'] = {}
			
			super_contigs[contig_number]['probes'][start] = {
				'stop': stop,
				'gene': gene
			}
			
	except(StopIteration):
		in_probe_file.close()

	return super_contigs


def generate_single_contig(
	in_fasta, 
	out_fasta,
	in_probe,
	out_probe,
	out_gene
):
	super_contigs = load_probes(in_probe)
	if not super_contigs:
		print (
			'!error: '
			'generate_single_contig(): '
			'cannot load probes from '+in_probe
		)
		return False

	#pprint(super_contigs)


	# load fasta file
	fasta_dict = SeqIO.index(in_fasta, 'fasta')
	#fasta = fasta_dict['supercont3.21']
	#print fasta.id
	#print fasta.seq

	# open new output fasta file
	try:
		out_fasta_file = open(out_fasta, 'w+')
	except IOError as err:
		print (
			'!error: '
			'generate_single_contig(): '
			'cannot open output fasta file: '+out_fasta
		)
		print err
		return False
	# write fasta header
	out_fasta_file.write('>chr1\n')

	# open new output probe file
	try:
		out_probe_file = open(out_probe, 'w+')
	except IOError as err:
		print (
			'!error: '
			'generate_single_contig(): '
			'cannot open output probe file: '+out_probe
		)
		print err
		return False
	# write probe header
	out_probe_file.write('chr\tstart\tstop\tname\n')


	# iterate through probes 
	cur_pos = 1
	gene_coords = {}

	for super_contig in sorted(super_contigs):
		fasta = fasta_dict[super_contigs[super_contig]['orig_name']]
		# print sequence to file
		out_fasta_file.write(str(fasta.seq))
		# translate probe coords and print to file
		for probe_start in sorted(super_contigs[super_contig]['probes']):
			probe = super_contigs[super_contig]['probes'][probe_start]

			start_pos = probe_start+cur_pos-1
			stop_pos = probe['stop']+cur_pos-1
			gene_name = probe['gene']

			out_probe_file.write(
				'1\t'+str(start_pos)+
				'\t'+str(stop_pos)+
				'\t'+gene_name+
				'\n'
			)
			if gene_name not in gene_coords:
				gene_coords[gene_name] = {}
				gene_coords[gene_name]['start'] = start_pos
				gene_coords[gene_name]['stop'] = stop_pos
			else:
				if start_pos < gene_coords[gene_name]['start']:
					gene_coords[gene_name]['start'] = start_pos
				if stop_pos > gene_coords[gene_name]['stop']:
					gene_coords[gene_name]['stop'] = stop_pos

		cur_pos += len(fasta.seq)

	fasta_dict.close()
	out_probe_file.close()
	out_fasta_file.close()

	# write gene coords
	try:
		out_gene_file = open(out_gene, 'w+')
	except IOError as err:
		print (
			'!error: '
			'generate_single_contig(): '
			'cannot open output gene coords file: '+out_gene
		)
		print err
		return False

	for gene in sorted(gene_coords):
		out_gene_file.write(
			gene+'\t'+
			str(gene_coords[gene]['start'])+'\t'+
			str(gene_coords[gene]['stop'])+'\n'
		)
	out_gene_file.close()
	

	return True


def main():
	(
		in_fasta,
		out_fasta,
		in_probe,
		out_probe,
		out_gene
	) = parse_args()

	if not generate_single_contig(in_fasta, out_fasta, in_probe, out_probe, out_gene):
		print (
			'!error: '
			'generate_single_contig() failed'
		)
		sys.exit(1)

	sys.exit(0)


if __name__ == "__main__":
	main()


