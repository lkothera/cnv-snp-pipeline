#!/bin/bash

# clean up generated reference files


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# remove files
rm -f ${BASE_REFERENCE}/${REFERENCE}-CDS.gtf
rm -f ${BASE_REFERENCE}/${REFERENCE}-combined.gtf
rm -f ${BASE_REFERENCE}/vectorbase_files/${VBFASTA}.fai
#rm -f ${BASE_REFERENCE}/gene_names_for_grep.txt
rm -f ${BASE_REFERENCE}/probes_final.txt
#rm -f ${BASE_REFERENCE}/probes_step1.txt
rm -f ${BASE_REFERENCE}/probes_full_step1.txt
#rm -f ${BASE_REFERENCE}/probes_step2.txt
rm -f ${BASE_REFERENCE}/probes_full_step2.txt
rm -f ${BASE_REFERENCE}/gene_coords.txt
rm -f ${BASE_REFERENCE}/bwa/supercontigs.fa*
rm -f ${BASE_REFERENCE}/bwa/supercontigs_combined.fa*
rm -f ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted.fa*
rm -f ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted.dict
rmdir ${BASE_REFERENCE}/bwa
rm -f ${BASE_REFERENCE}/supercontigs_uniq.txt
