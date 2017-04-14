#!/bin/bash


# general variables
BASE=/scicomp/home/ktr2/Projects/Mosquito2
REFERENCE=aedes2
SAMPLES=aedes2

# variables for reference generation
VBGTF=aedes-aegypti-liverpoolbasefeaturesaaegl33.gtf
VBFASTA=aedes-aegypti-liverpoolscaffoldsaaegl3.fa
BASE_REFERENCE=${BASE}/reference/${REFERENCE}

# variables for sample QC processing
FASTQ_DIR=fastq_orig
FAQCS_DIR=fastq_faqcs
BASE_SAMPLES=${BASE}/samples/${SAMPLES}

# variables for conifer
RUN_NAME=run-aedes
BWA_REFERENCE=${BASE_REFERENCE}/bwa/supercontigs_combined_formatted.fa

# variables for feature selection
COMPARE=ResistancePhenos_QonlyAsOf11-8-16.txt
