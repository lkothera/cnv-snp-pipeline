#!/usr/bin/env python

# -----
# Title: assoc_variants.py
# Author: John Phan (John.Phan@csra.com)
# Last Modified: 2017/07/27
# -----

# associate variants in vcf to groups


import os
import re
import glob
import sys
import argparse
from pprint import pprint
from subprocess import (PIPE, Popen)
import vcf
import scipy.stats as stats
import statsmodels.api as sm 
from Bio import SeqIO
from sortedcontainers import SortedDict
from string import maketrans

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
	parser.add_argument(
		'--fasta',
		type=str,
		required=True,
		help='reference fasta file'
	)
	parser.add_argument(
		'--gtf',
		type=str,
		required=True,
		help='probe gtf file'
	)
	args = parser.parse_args()
	return (
		args.vcf,
		args.compare,
		args.column,
		args.fasta,
		args.gtf
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

def load_gtf(
	gtf_path
):

	# open gtf file
	try:
		gtf_file = open(gtf_path, 'rU')
	except IOError as err:
		print (
			'!error: '
			'load_gtf(): '
			'cannot open gtf file: '+gtf_path
		)
		print err
		return False

	# initialize dict data structure
	gtf = {}
	genes = {}
	try:
		while True:
			line = gtf_file.next().strip()
			line_elems = re.split(r'\t', line)
			contig = line_elems[0]
			start = int(line_elems[3])
			stop = int(line_elems[4])
			strand = line_elems[6]
			frame = int(line_elems[7])
			properties = line_elems[8]

			regex = re.compile(r"gene_id \"(.*?)\"; transcript_id \"(.*?)\"; exon_number \"(.*?)\";")
			matches = regex.search(properties).groups()
			gene_id = matches[0]
			transcript_id = matches[1]
			exon_number = matches[2]

			if contig not in gtf:
				gtf[contig] = SortedDict()

			gtf[contig][start] = {
				'stop': stop,
				'gene_id': gene_id,
				'transcript_id': transcript_id,
				'exon_number': exon_number,
				'frame': frame
			}

			# keep track of gene start/stop
			if gene_id not in genes:
				genes[gene_id] = {
					'start': start,
					'stop': stop,
					'contig': contig,
					'strand': strand
				}
			else:
				if start < genes[gene_id]['start']:
					genes[gene_id]['start'] = start
				if stop > genes[gene_id]['stop']:
					genes[gene_id]['stop'] = stop

	except(StopIteration):
		gtf_file.close()

	# create sorted dict for gene start/stop
	genes_sorted = {}
	for gene_id in genes:
		contig = genes[gene_id]['contig']
		if contig not in genes_sorted:
			genes_sorted[contig] = SortedDict()
		genes_sorted[contig][genes[gene_id]['start']] = {
			'gene_id': gene_id,
			'stop': genes[gene_id]['stop']
		}

	return (gtf, genes_sorted, genes)


def get_variant_info(
	gtf_dict,
	genes_sorted,
	genes,
	fasta_dict,
	vcf_record
):

	var_pos = vcf_record.POS
	gtf_key = 0
	gene_id = ''

	# determine if variant position corresponds to a gene, intron, or integenic region
	index = gtf_dict[vcf_record.CHROM].bisect_left(var_pos)
	if index == 0:
		# variant position is before, or exactly on an exon start location
		if var_pos in gtf_dict[vcf_record.CHROM]:
			# variant is exactly on start location
			gtf_key = gtf_dict[vcf_record.CHROM].iloc[index]
			gene_id = gtf_dict[vcf_record.CHROM][gtf_key]['gene_id']
		else:
			# variant is before, non-cds
			gene_id = 'Non-CDS'
	else:
		# variant position is after exon location
		gtf_key = gtf_dict[vcf_record.CHROM].iloc[index-1]
		if var_pos <= gtf_dict[vcf_record.CHROM][gtf_key]['stop']:
			# if inside of the exon's region
			gene_id = gtf_dict[vcf_record.CHROM][gtf_key]['gene_id']
		else:
			# outside of exon's region
			gene_id = 'Non-CDS'

	if gene_id == 'Non-CDS':
		# determine if in an intron, or intergenic region
		index = genes_sorted[vcf_record.CHROM].bisect_left(var_pos)
		if index > 0:
			gene_key = genes_sorted[vcf_record.CHROM].iloc[index-1]
			if var_pos <= genes_sorted[vcf_record.CHROM][gene_key]['stop']:
				# append intron gene id
				gene_id += ' (Intron '+genes_sorted[vcf_record.CHROM][gene_key]['gene_id']+')'

	

	codon_seq = 'N/A'
	alt_codon_seq = 'N/A'
	codon_pos = 0
	pattern = re.compile("^Non-CDS")
	if not pattern.match(gene_id):
		# only if coding sequence
		
		# extra fasta sequence for entire "chromosome"
		fasta_seq = fasta_dict[vcf_record.CHROM]

		if genes[gene_id]['strand'] == '+':
			# if forward strand, calculate codon info from the beginning

			# beginning of first whole codon for this exon
			cds_start = gtf_key+gtf_dict[vcf_record.CHROM][gtf_key]['frame']
			# beginning of codon that contains the variant position
			codon_start = (var_pos-cds_start)//3*3+cds_start
			# position of the variant in the codon (1-based)
			codon_pos = (var_pos-cds_start)%3+1
			# extract codon sequence
			codon_seq = str(fasta_seq.seq[(codon_start-1):(codon_start+2)].upper())
			# derive alt. codon
			if len(vcf_record.REF) > 1:
				# all alternates are indels
				alt_codon_seq = '[indel]'
			else:
				alt_codon_seq = ','.join([ codon_seq[:codon_pos-1]+str(x)+codon_seq[codon_pos-1+len(str(x)):] if len(str(x)) == 1 else '[indel]' for x in vcf_record.ALT ])
		else:
			# if reverse strand, calculate codon from the end

			# beginning of first whole codon for this exon (counted from the end)
			cds_start = gtf_dict[vcf_record.CHROM][gtf_key]['stop']-gtf_dict[vcf_record.CHROM][gtf_key]['frame']
			# beginning of codon that contains the variant position
			codon_start = cds_start-(cds_start-var_pos)//3*3
			# position of the variant in the codon (1-based)
			codon_pos = (cds_start-var_pos)%3+1
			# extract codon sequence
			codon_seq = str(fasta_seq.seq[(codon_start-3):(codon_start)].upper())
			# derive alt. codon
			if len(vcf_record.REF) > 1:
				# all alternates are indels
				alt_codon_seq = '[indel]'
			else:
				alt_codon_seq = ','.join([ codon_seq[:(3-codon_pos)]+str(x)+codon_seq[(3-codon_pos)+len(str(x)):] if len(str(x)) == 1 else '[indel]' for x in vcf_record.ALT ])

			# reverse complement the codon sequences if not indels
			codon_seq = codon_seq + ' ('+codon_seq[::-1].translate(maketrans('ATGC','TACG'))+')'
			alt_codon_seq = (
				alt_codon_seq + 
				' (' + 
					','.join(
						[ x[::-1].translate(maketrans('ATGC','TACG')) if x != '[indel]' else '[indel]' for x in alt_codon_seq.split(',') ]
					) +
				')'
			)

	return (
		gene_id,
		codon_seq,
		alt_codon_seq,
		codon_pos
	) 

def main():
	(
		vcf_path,
		compare,
		column,
		fasta,
		gtf
	) = parse_args()

	compare_dict = read_compare(compare, column)
	if not compare_dict:
		print (
			"!error: "
			"assoc_variants.py: "
			"read_compare() failed"
		)
		sys.exit(1)


	# load gtf
	(gtf_dict, genes_sorted, genes) = load_gtf(gtf)

	
	# load fasta
	fasta_dict = SeqIO.index(fasta, 'fasta')

	print "ID\tPOS\tGene\tREF\tALT\tCodon\tAlt. Codon\tCodon Pos.\t# Called\t# Heterozygous\t# Homozygous Alt. Allele\t# Homozygous Ref. Allele\t# Ref. Susceptible (Group 1)\t# Ref. Resistant (Group 2)\t# Heterozygous Alt. Susceptible (Group 1)\t# Homozygous Alt. Susceptible (Group 1)\t# Heterozygous Alt. Resistant (Group 2)\t# Homozygous Alt. Resistant (Group 2)\tP-Val., Fisher Exact (Method 1)\tFE Statistic (Method 1)\tOdds Ratio (Method 1)\tSensitivity (Method 1)\tSpecificity (Method 1)\tAUC (Average Sens,Spec; Method 1)\tAccuracy (Method 1)\tP-Val., Fisher Exact (Method 2)\tFE Statistic (Method 2)\tOdds Ratio (Method 2)\tSensitivity (Method 2)\tSpecificity (Method 2)\tAUC (Average Sens,Spec; Method 2)\tAccuracy (Method 2)\tP-Val., Cochran-Armitage\tCochran-Armitage Statistic\tHeterozygous Alt. Susceptible Samples\tHomozygous Alt. Susceptible Samples\tHeterozygous Alt. Resistant Samples\tHomozygous Alt. Resistant Samples"

	vcf_reader = vcf.Reader(open(vcf_path, 'r'));
	for record in vcf_reader:

		#print record.REF
		#print record.CHROM
		res_ref = 0
		res_alt_hom = 0
		res_alt_het = 0
		sus_ref = 0
		sus_alt_hom = 0
		sus_alt_het = 0

		num_het = 0
		num_hom_ref = 0
		num_hom_alt = 0
		sus_alt_hom_list = []
		sus_alt_het_list = []
		res_alt_hom_list = []
		res_alt_het_list = []

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
					# alt type
					if compare_dict[sample_name] == 1:
						# susceptible
						sus_alt_het += 1
						sus_alt_het_list.append(sample_name)
					else:
						res_alt_het += 1
						res_alt_het_list.append(sample_name)
				else:
					num_hom_alt += 1
					# alt type
					if compare_dict[sample_name] == 1:
						# susceptible
						sus_alt_hom += 1
						sus_alt_hom_list.append(sample_name)
					else:
						res_alt_hom += 1
						res_alt_hom_list.append(sample_name)

		# calculate FE based on grouping count of hom alt and het alt
		oddsratio, pvalue = stats.fisher_exact([[res_ref, sus_ref], [res_alt_hom+res_alt_het, sus_alt_hom+sus_alt_het]])
		stat = 0
		try:
			stat = res_ref/float(res_ref+sus_ref)-(res_alt_hom+res_alt_het)/float(res_alt_hom+res_alt_het+sus_alt_hom+sus_alt_het)
		except ZeroDivisionError:
			# p-value will be non-significant anyway
			stat = 0

		# calculate FE based on individual counts of allele copies in hom alt (2) and het alt (1)
		oddsratio2, pvalue2 = stats.fisher_exact([[2*res_ref, 2*sus_ref], [2*res_alt_hom+res_alt_het, 2*sus_alt_hom+sus_alt_het]])
		stat2 = 0
		try:
			stat2 = 2*res_ref/float(2*res_ref+2*sus_ref)-(2*res_alt_hom+res_alt_het)/float(2*res_alt_hom+res_alt_het+2*sus_alt_hom+sus_alt_het)
		except ZeroDivisionError:
			stat2 = 0

		# calculate cochran-armitage test
		catt_table = sm.stats.Table([[res_ref,sus_ref], [res_alt_het,sus_alt_het], [res_alt_hom,sus_alt_hom]])
		catt_result = catt_table.test_ordinal_association()

		# sensitivity & specificity using presence/absense of alt alleles as predictor
		sens1 = 0
		spec1 = 0
		auc1 = 0
		acc1 = 0

		sens2 = 0
		spec2 = 0
		auc2 = 0
		acc2 = 0
		
		#if float(pvalue) < thresh: # check for method 1 of FE statistic
		if stat < 0: # allele presence indicates resistance
			try:
				sens1 = (res_alt_hom+res_alt_het)/float(res_alt_hom+res_alt_het+res_ref)
			except ZeroDivisionError:
				pass
			try:
				spec1 = (sus_ref)/float(sus_ref+sus_alt_hom+sus_alt_het)
			except ZeroDivisionError:
				pass
			auc1 = (sens1+spec1)/float(2)
			try:
				acc1 = (res_alt_hom+res_alt_het+sus_ref)/float(num_het+num_hom_ref+num_hom_alt)
			except ZeroDivisionError:
				pass
		else: # allele presence indicates susceptible
			try:
				sens1 = (res_ref)/float(res_ref+res_alt_hom+res_alt_het)
			except ZeroDivisionError:
				pass
			try:
				spec1 = (sus_alt_hom+sus_alt_het)/float(sus_alt_hom+sus_alt_het+sus_ref)
			except ZeroDivisionError:
				pass
			auc1 = (sens1+spec1)/float(2)
			try:
				acc1 = (sus_alt_hom+sus_alt_het+res_ref)/float(num_het+num_hom_ref+num_hom_alt)
			except ZeroDivisionError:
				pass
		#if float(pvalue2) < thresh: # check for method 2 of FE statistic
		if stat2 < 0: # allele presence indicates resistance
			try:
				sens2 = (res_alt_hom+res_alt_het)/float(res_alt_hom+res_alt_het+res_ref)
			except ZeroDivisionError:
				pass
			try:
				spec2 = (sus_ref)/float(sus_ref+sus_alt_hom+sus_alt_het)
			except ZeroDivisionError:
				pass
			auc2 = (sens2+spec2)/float(2)
			try:
				acc2 = (res_alt_hom+res_alt_het+sus_ref)/float(num_het+num_hom_ref+num_hom_alt)
			except ZeroDivisionError:
				pass
		else: # allele presence indicates susceptible
			try:
				sens2 = (res_ref)/float(res_ref+res_alt_hom+res_alt_het)
			except ZeroDivisionError:
				pass
			try:
				spec2 = (sus_alt_hom+sus_alt_het)/float(sus_alt_hom+sus_alt_het+sus_ref)
			except ZeroDivisionError:
				pass
			auc2 = (sens2+spec2)/float(2)
			try:
				acc2 = (sus_alt_hom+sus_alt_het+res_ref)/float(num_het+num_hom_ref+num_hom_alt)
			except ZeroDivisionError:
				pass


		# lookup corresponding gene and codon info
		(
			gene_id,
			codon_seq,
			alt_codon_seq,
			codon_pos
		) = get_variant_info(gtf_dict, genes_sorted, genes, fasta_dict, record)


		print (
			record.CHROM + "\t" +
			str(record.POS) + "\t" +
			gene_id + "\t" +
			record.REF + "\t" +
			",".join(map(str, record.ALT)) + "\t" +
			codon_seq + "\t" +
			alt_codon_seq + "\t" +
			str(codon_pos) + "\t" +
			str(num_het+num_hom_ref+num_hom_alt) + "\t" +
			str(num_het) + "\t" +
			str(num_hom_alt) + "\t" +
			str(num_hom_ref) + "\t" +
			str(sus_ref) + "\t" +
			str(res_ref) + "\t" +
			str(sus_alt_het) + "\t" +
			str(sus_alt_hom) + "\t" +
			str(res_alt_het) + "\t" +
			str(res_alt_hom) + "\t" +
			str(pvalue) + "\t" +
			str(stat) + "\t" +
			str(oddsratio) + "\t" +
			str(sens1) + "\t" +
			str(spec1) + "\t" +
			str(auc1) + "\t" +
			str(acc1) + "\t" +
			str(pvalue2) + "\t" +
			str(stat2) + "\t" +
			str(oddsratio2) + "\t" +
			str(sens2) + "\t" +
			str(spec2) + "\t" +
			str(auc2) + "\t" +
			str(acc2) + "\t" +
			str(catt_result.pvalue) + "\t" +
			str(catt_result.zscore) + "\t" +
			",".join(sus_alt_het_list) + "\t" +
			",".join(sus_alt_hom_list) + "\t" +
			",".join(res_alt_het_list) + "\t" +
			",".join(res_alt_hom_list)
		)


if __name__ == "__main__":
	main()


