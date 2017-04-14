#!/bin/bash

# general variables
BASE=/scicomp/home/ktr2/Projects/Mosquito3
REFERENCE=culex
SAMPLES=culex

# variables for reference generation
VBGTF=culex-quinquefasciatus-johannesburgbasefeaturescpipj22.gtf
VBFASTA=culex-quinquefasciatus-johannesburgscaffoldscpipj2.fa
BASE_REFERENCE=${BASE}/reference/${REFERENCE}

# variables for sample QC processing
FASTQ_DIR=fastq_orig
FAQCS_DIR=fastq_faqcs
BASE_SAMPLES=${BASE}/samples/${SAMPLES}

# variables for conifer
RUN_NAME=run-culex
BWA_REFERENCE=${BASE_REFERENCE}/bwa/supercontigs_combined_formatted.fa

# variables for feature selection
COMPARE=ResistancePhenos_QonlyAsOf11-8-16.txt
