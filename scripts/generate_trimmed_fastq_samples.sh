#!/bin/bash

# QC analysis and trimming of fastq samples


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# load dependencies
source /etc/profile.d/modules.sh
module load FaQCs_OAMD/1.34

# create subdirectories
mkdir -p ${BASE_SAMPLES}/${FAQCS_DIR}
mkdir -p ${BASE_SAMPLES}/${FAQCS_DIR}_only

# create list of all samples
find ${BASE_SAMPLES}/${FASTQ_DIR} | grep '\.fastq' | sed 's/.*\/\(.*\)\.fastq/\1/' | sort > ${BASE_SAMPLES}/sample_names.txt

# run FaQCs.pl on each sample
ls ${BASE_SAMPLES}/${FASTQ_DIR} | xargs --verbose -n 1 -I % sh -c "FaQCs.pl -u ${BASE_SAMPLES}/${FASTQ_DIR}/% -q 30 -d ${BASE_SAMPLES}/${FAQCS_DIR}/%"

# copy trimmed fastq files
find ${BASE_SAMPLES}/${FAQCS_DIR} | grep fastq.gz | sed 's/.*\/\(.*\)\.fastq\/.*/\1/' | xargs --verbose -n 1 -I % sh -c "cp ${BASE_SAMPLES}/${FAQCS_DIR}/%.fastq/QC.unpaired.trimmed.fastq.gz ${BASE_SAMPLES}/${FAQCS_DIR}_only/%.fastq.gz"

# unzip trimmed fastq files
gunzip ${BASE_SAMPLES}/${FAQCS_DIR}_only/*

