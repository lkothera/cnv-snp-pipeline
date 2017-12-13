#!/usr/bin/env python

# -----
# Title: combine_rpkms.py
# Author: John Phan (John.Phan@csra.com, ktr2@cdc.gov)
# -----

# combine rpkm files from conifer, use combined exports file as guide for 
# inclusion of probes


import os
import glob
import re
import sys
import argparse
from pprint import pprint
from subprocess import (PIPE, Popen)
import math


def parse_args():
    # parse command line arguments and return
    # as list of argument values
    parser = argparse.ArgumentParser(
        description='combine_rpkms.py'
    )
    parser.add_argument(
        '--in_dir',
        type=str,
        required=True,
        help='input directory'
    )
    parser.add_argument(
        '--combined_exports_file',
        type=str,
        required=True,
        help='combined exports file'
    )
    parser.add_argument(
        '--probes_file',
        type=str,
        required=True,
        help='probes file'
    )
    parser.add_argument(
        '--out',
        type=str,
        required=True,
        help='output data file'
    )
    parser.add_argument(
        '--out_gene',
        type=str,
        required=True,
        help='output data file for combined gene info'
    )
    parser.add_argument(
        '--log',
        type=str,
        required=True,
        help='log2 values?'        
    )
    args = parser.parse_args()
    return (
        args.in_dir,
        args.combined_exports_file,
        args.probes_file,
        args.out,
        args.out_gene,
        args.log
    )


def read_probes_file(
    probes_file
):

    # open probes file
    try:
        probes_f = open(probes_file, 'rU')
    except IOError as err:
        print(
            '!error: read_probes_file(): '
            'cannot open probes file: '+probes_file
        )
        print(err)
        return False

    probes = {}

    # skip header
    line = probes_f.next()
    # read remaining line by line
    line_num = 1
    try:
        while True:
            line = probes_f.next().strip()
            line_elems = re.split(r'\t', line)
            probe_start = line_elems[1]
            gene_name = line_elems[3]

            if gene_name not in probes:
                probes[gene_name] = {}

            # store the line (or probe) number for each probe
            probes[gene_name][probe_start] = line_num
            line_num += 1

    except(StopIteration):
        probes_f.close()

    return probes
    
def read_export_probes(
    exports_file
):
    # open export file
    try:
        exports_f = open(exports_file, 'rU')
    except IEError as err:
        print(
            '!error read_export_probes(): '
            'cannot open export file: '+exports_file
        )
        print(err)
        return False

    exports = []

    # skip header
    line = exports_f.next()
    # read remaining line by line
    try:
        while True:
            line = exports_f.next().strip()
            line_elems = re.split(r'\t', line)
            probe = line_elems[0]
            probe_elems = re.split(r'_', probe)

            exports.append({
                'gene': probe_elems[0],
                'start': probe_elems[1]
            })
    except(StopIteration):
        exports_f.close()

    return exports


def read_rpkm_files(
    in_dir
):

    rpkm_data = {}

    # get list of fastq files
    dirlist = glob.glob(in_dir+'/*.rpkm.txt')

    for rpkm_file_name in dirlist:
        basename = os.path.basename(rpkm_file_name)
        sample_name = os.path.splitext(basename)[0]
        sample_name = os.path.splitext(sample_name)[0]

        # open bed file
        try:
            rpkm_file = open(rpkm_file_name, 'rU')
        except IOError as err:
            print (
                '!error: '
                'read_rpkm_files(): '
                'cannot open rpkm file: '+rpkm_file_name
            )
            print err
            return False


        try:
            while True:
                line = rpkm_file.next().strip()
                line_elems = re.split(r'\t', line)

                probe_num = line_elems[0]
                rpkm_val = line_elems[2]

                if sample_name not in rpkm_data:
                    rpkm_data[sample_name] = {}

                rpkm_data[sample_name][probe_num] = rpkm_val

        except(StopIteration):
            rpkm_file.close()

    return rpkm_data


def main():
    (
        in_dir,
        combined_exports_file,
        probes_file,
        out,
        out_gene,
        log
    ) = parse_args()

    # read probes file
    probes = read_probes_file(probes_file)
    #pprint(probes)


    # read export file probe names
    export_probes = read_export_probes(combined_exports_file)
    #pprint(export_probes)


    # read rpkm files
    data = read_rpkm_files(in_dir)
    #pprint(data)


    try:
        out_file = open(out, 'w+')
    except IOError as err:
        print (
            '!error: '
            'combine_rpkms.py: '
            'cannot open output file: '+out
        )
        print(err)
        sys.exit(1)

    # print header
    for sample in sorted(data):
        out_file.write('\t'+sample)
    out_file.write('\n')

    # write all probe info
    for probe in sorted(probes):
        for coord in sorted(probes[probe]):
            out_file.write(probe+'_'+coord)
            line_num = str(probes[probe][coord])
            for sample in sorted(data):
                value = float(data[sample][line_num])
                if log=='true':
                    value = math.log(value+1, 2)
                out_file.write('\t'+str(value))
            out_file.write('\n')

    out_file.close()

    # write gene info
    try:
        out_file = open(out_gene, 'w+')
    except IOError as err:
        print(
            '!error: combine_rpkms.py: '
            'cannot open output file: '+out_gene
        )
        print(err)
        sys.exit(1)

    # print header
    for sample in sorted(data):
        out_file.write('\t'+sample)
    out_file.write('\n')

    # write all probe info
    for probe in sorted(probes):
        out_file.write(probe)
        probe_data = {}
        for coord in sorted(probes[probe]):
            line_num = str(probes[probe][coord])
            for sample in sorted(data):
                if sample not in probe_data:
                    probe_data[sample] = []
                rpkm_value = float(data[sample][line_num])
                if rpkm_value > 0:
                    probe_data[sample].append(rpkm_value)
        for sample in sorted(probe_data):
            if len(probe_data[sample]) > 0:
                value = float(sum(probe_data[sample])/len(probe_data[sample]))
                if log=='true':
                    value = math.log(value+1, 2)
                out_file.write('\t'+str(value))
            else:
                out_file.write('\t0.0')
        out_file.write('\n')

    out_file.close()

    sys.exit(0)


if __name__ == "__main__":
    main()


