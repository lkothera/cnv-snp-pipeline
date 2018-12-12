#!/bin/bash

# clean up generated reference files


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# remove files
rm -f ${BASE_REFERENCE}/${REFERENCE}-CDS.gtf
rm -f ${BASE_REFERENCE}/${REFERENCE}-exons.gtf
rm -f ${BASE_REFERENCE}/${REFERENCE}-combined_gatk.gtf
rm -f ${BASE_REFERENCE}/vectorbase_files/${VBFASTA}.fai
rm -f ${BASE_REFERENCE}/gene_names_for_grep.txt
rm -f ${BASE_REFERENCE}/probes_final_conifer.txt
rm -f ${BASE_REFERENCE}/probes_final_gatk.txt
rm -f ${BASE_REFERENCE}/probes_step1_conifer.txt
rm -f ${BASE_REFERENCE}/probes_step1_gatk.txt
rm -f ${BASE_REFERENCE}/probes_step2_conifer.txt
rm -f ${BASE_REFERENCE}/probes_step2_gatk.txt
rm -f ${BASE_REFERENCE}/gene_coords_conifer.txt
rm -f ${BASE_REFERENCE}/gene_coords_gatk.txt
rm -f ${BASE_REFERENCE}/bwa/*
rmdir ${BASE_REFERENCE}/bwa
rm -f ${BASE_REFERENCE}/supercontigs_uniq_conifer.txt
rm -f ${BASE_REFERENCE}/supercontigs_uniq_gatk.txt
