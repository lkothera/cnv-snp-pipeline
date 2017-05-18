#!/bin/bash

# clean up generated reference files


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# remove files
rm ${BASE_REFERENCE}/${REFERENCE}-CDS.gtf
rm ${BASE_REFERENCE}/vectorbase_files/${VBFASTA}.fai
rm ${BASE_REFERENCE}/gene_names_for_grep.txt
rm ${BASE_REFERENCE}/probes_final.txt
rm ${BASE_REFERENCE}/probes_step1.txt
rm ${BASE_REFERENCE}/probes_step2.txt
rm ${BASE_REFERENCE}/gene_coords.txt
rm ${BASE_REFERENCE}/bwa/supercontigs.fa*
rm ${BASE_REFERENCE}/bwa/supercontigs_combined.fa*
rm ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted.fa*
rm ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted.dict
rmdir ${BASE_REFERENCE}/bwa
rm ${BASE_REFERENCE}/supercontigs_uniq.txt
