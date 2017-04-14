#!/bin/bash


# run ion-torrent GATK pipeline


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# load dependencies
source /etc/profile.d/modules.sh
module load Python/2.7
module load perl/5.16.1-MT
module load picard/1.128
module load gatk/3.6
module load samtools/1.3.1
module load tmap/3.4.1
module load bcftools/1.2
module load vcftools/0.1.14


# create run directory 
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/gatk

# map fastq files to reference
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/gatk/tmap
${BASE}/scripts/tmap_align.py --in_dir ${BASE}/samples/${SAMPLES}/${FAQCS_DIR}_only --out_dir ${BASE}/pipeline-runs/${RUN_NAME}/gatk/tmap --base_dir ${BASE} --ext fastq --database ${BWA_REFERENCE} 

# run pre_variant script
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/gatk/pre_variant
${BASE}/scripts/pre_variant.py --in_dir ${BASE}/pipeline-runs/${RUN_NAME}/gatk/tmap --out_dir ${BASE}/pipeline-runs/${RUN_NAME}/gatk/pre_variant --base_dir ${BASE}

# realign indels
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/gatk/indel_realign
${BASE}/scripts/indel_realign.py --in_dir ${BASE}/pipeline-runs/${RUN_NAME}/gatk/pre_variant --out_dir ${BASE}/pipeline-runs/${RUN_NAME}/gatk/indel_realign --reference ${BWA_REFERENCE}

# create bam list
find ${BASE}/pipeline-runs/${RUN_NAME}/gatk/indel_realign | grep rg.bam$ > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/bam.list

# gatk haplotype caller
GenomeAnalysisTK.sh -T HaplotypeCaller -l INFO -I ${BASE}/pipeline-runs/${RUN_NAME}/gatk/bam.list -stand_call_conf 30 -stand_emit_conf 10 -R ${BWA_REFERENCE} -o ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.vcf -nct 8

# extract header from intermediate vcf results
awk '{if(($1~/#/)) print $0}' ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.vcf > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.head

# mpileup of all bams
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant
samtools mpileup -f ${BWA_REFERENCE} `find ${BASE}/pipeline-runs/${RUN_NAME}/gatk/indel_realign | grep rg.bam$` > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.mpileup

# filter mpileup to remove low depth bases
cat ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.mpileup | awk '{if($4 >= 10) print $0}' > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants_filtered.mpileup

# load a different perl - pre-variant perl script only works w/ MT perl while post-variant perl scripts only work w/ newer perl
module unload perl/5.16.1-MT
module load perl/5.22.1

# extract variants
${BASE}/scripts/pileup_analyse.pl ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants_filtered.mpileup ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.snp.indel

# filter biases
${BASE}/scripts/filter_indel_bias.pl -vcf ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.vcf -pileup ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.snp.indel -o ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.indel1

# filter homopolymer indels
${BASE}/scripts/filter_indel_homo2.pl -vcf ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.indel1 -g ${BWA_REFERENCE} -o ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.final.indel

# create final indel vcf
cat ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.head ${BASE}/pipeline-runs/${RUN_NAME}/gatk/post_variant/variants.final.indel > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.indel.vcf

# extract snps into separate vcf
bcftools view --types snps ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.vcf > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.snp.vcf

# format into easy-to-read tables
vcf-to-tab < ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.snp.vcf > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.snp_table.txt
vcf-to-tab < ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.indel.vcf > ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.indel_table.txt


