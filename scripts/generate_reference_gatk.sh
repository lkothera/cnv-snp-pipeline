#!/bin/bash

# generate a reference and probe file for GATK using VectorBase source files

# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# load dependencies
source /etc/profile.d/modules.sh
module load bwa/0.7.12
module load samtools/1.3.1
module load Python/2.7.12
module load seqtk/1.0
module load picard/1.128
module load tmap/3.4.1

# extract only coding sequences from the gtf file
grep "\sCDS\s\|\sstop_codon\s" ${BASE_REFERENCE}/vectorbase_files/${VBGTF} > ${BASE_REFERENCE}/${REFERENCE}-CDS.gtf

# get list of unique supercontigs
cat ${BASE_REFERENCE}/${REFERENCE}-CDS.gtf | awk '{print $1}' | sort | uniq > ${BASE_REFERENCE}/supercontigs_uniq_gatk.txt

# index fasta file
samtools faidx ${BASE_REFERENCE}/vectorbase_files/${VBFASTA}

# generate supercontigs file
mkdir -p ${BASE_REFERENCE}/bwa
cat ${BASE_REFERENCE}/supercontigs_uniq_gatk.txt | xargs --verbose -n 1 -I % sh -c "samtools faidx ${BASE_REFERENCE}/vectorbase_files/${VBFASTA} % >> ${BASE_REFERENCE}/bwa/supercontigs_gatk.fa"

# generate probes file
cat ${BASE_REFERENCE}/${REFERENCE}-CDS.gtf | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$12"\t"$14}' > ${BASE_REFERENCE}/probes_step1_gatk.txt
cat ${BASE_REFERENCE}/probes_step1_gatk.txt | perl -pe 's/"(.*?)";/\1/g' > ${BASE_REFERENCE}/probes_step2_gatk.txt

# combine supercontigs and remap probe coordinates
${BASE}/scripts/generate_reference_contigs.py --in_fasta ${BASE_REFERENCE}/bwa/supercontigs_gatk.fa --out_fasta ${BASE_REFERENCE}/bwa/supercontigs_combined_gatk.fa --in_probe ${BASE_REFERENCE}/probes_step2_gatk.txt --out_probe ${BASE_REFERENCE}/probes_final_gatk.txt --out_gene ${BASE_REFERENCE}/gene_coords_gatk.txt --out_gtf ${BASE_REFERENCE}/${REFERENCE}-combined_gatk.gtf

# reformat fasta file
seqtk seq -l 100 ${BASE_REFERENCE}/bwa/supercontigs_combined_gatk.fa > ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.fa

# create bwa index
bwa index ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.fa

# create tmap index
tmap index -f ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.fa

# index fasta file
samtools faidx ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.fa

# create picard sequence dictionary for gatk
picard CreateSequenceDictionary REFERENCE=${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.fa OUTPUT=${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.dict
