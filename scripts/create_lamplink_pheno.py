#!/usr/bin/env python

# -----
# Title: create_lamplink_pheno.py
# Author: John Phan (John.Phan@csra.com)
# -----

# create a phenotype file for lamplink


import os
import re
import sys
import argparse
from pprint import pprint


def parse_args():

    # parse command line arguments and return
    # as list of argument values
    parser = argparse.ArgumentParser(
        description='create_lamplink_pheno.py'
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


def print_pheno_data(
    sample_to_label_map
):

    # print header
    sys.stdout.write('ID\tFID\tphe\n')

    # print each sample
    for sample in sample_to_label_map:
        sys.stdout.write('%s.rmDup_sorted\t%s.rmDup_sorted\t%s\n' % (
            sample, sample, str(sample_to_label_map[sample])
        ))
    
    return True


def main():

    (
        label_path,
        column
    ) = parse_args()

    sample_to_label_map = load_sample_to_label_map(label_path, column)
    if not sample_to_label_map:
        print (
            '!error: '
            'create_lamplink_pheno.py: '
            'cannot load sample to label map'
        )
        sys.exit(1)

    print_pheno_data(sample_to_label_map)

    sys.exit(0)


if __name__ == "__main__":
    main()

