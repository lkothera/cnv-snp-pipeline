#!/bin/bash

# generate a reference and probe file for conifer using VectorBase source files


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

# extract only exons from the gtf file
grep "\sexon\s" ${BASE_REFERENCE}/vectorbase_files/${VBGTF} > ${BASE_REFERENCE}/${REFERENCE}-exons.gtf

# create string for grep OR
cat ${BASE_REFERENCE}/gene_names.txt | paste -sd"," | sed 's/,/\\|/g' > ${BASE_REFERENCE}/gene_names_for_grep.txt

# get list of unique supercontigs
grep `cat ${BASE_REFERENCE}/gene_names_for_grep.txt` ${BASE_REFERENCE}/${REFERENCE}-exons.gtf | awk '{print $1}' | sort | uniq > ${BASE_REFERENCE}/supercontigs_uniq_conifer.txt

# index fasta file
samtools faidx ${BASE_REFERENCE}/vectorbase_files/${VBFASTA}

# generate supercontigs file
mkdir -p ${BASE_REFERENCE}/bwa
cat ${BASE_REFERENCE}/supercontigs_uniq_conifer.txt | xargs --verbose -n 1 -I % sh -c "samtools faidx ${BASE_REFERENCE}/vectorbase_files/${VBFASTA} % >> ${BASE_REFERENCE}/bwa/supercontigs_conifer.fa"

# generate probes file for conifer
grep `cat ${BASE_REFERENCE}/gene_names_for_grep.txt ` ${BASE_REFERENCE}/${REFERENCE}-exons.gtf | awk '{print $1"\t"$4"\t"$5"\t"$10}' > ${BASE_REFERENCE}/probes_step1_conifer.txt
echo "chr	start	stop	name" > ${BASE_REFERENCE}/probes_step2_conifer.txt
cat ${BASE_REFERENCE}/probes_step1_conifer.txt | sed 's/\"\(.*\)\";/\1/' >> ${BASE_REFERENCE}/probes_step2_conifer.txt

# combine supercontigs and remap probe coordinates
${BASE}/scripts/generate_single_contig.py --in_fasta ${BASE_REFERENCE}/bwa/supercontigs_conifer.fa --out_fasta ${BASE_REFERENCE}/bwa/supercontigs_combined_conifer.fa --in_probe ${BASE_REFERENCE}/probes_step2_conifer.txt --out_probe ${BASE_REFERENCE}/probes_final_conifer.txt --out_gene ${BASE_REFERENCE}/gene_coords_conifer.txt

# reformat fasta file
seqtk seq -l 100 ${BASE_REFERENCE}/bwa/supercontigs_combined_conifer.fa > ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_conifer.fa

# create bwa index
bwa index ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_conifer.fa

# create tmap index
tmap index -f ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_conifer.fa

# index fasta file
samtools faidx ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_conifer.fa

# create picard sequence dictionary
picard CreateSequenceDictionary REFERENCE=${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_conifer.fa OUTPUT=${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_conifer.dict
