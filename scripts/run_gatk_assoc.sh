#!/bin/bash

# associate SNPs in vcf file with groups

# parse arguments
for i in "$@"
do
case $i in
    -n=*|--name=*)
    NAME="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--column=*)
    COLUMN="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

echo "NAME = ${NAME}"
echo "COLUMN = ${COLUMN}"

# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# load dependencies
source /etc/profile.d/modules.sh
module load Python/2.7

# create features folder
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/gatk_features

# run feature assoc script
#${BASE}/scripts/assoc_variants2.py --vcf ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.snp.vcf --compare ${BASE}/samples/${SAMPLES}/${COMPARE} --column ${COLUMN} > ${BASE}/pipeline-runs/${RUN_NAME}/gatk_features/variants.snp_assoc_${NAME}.txt
${BASE}/scripts/assoc_variants.py --vcf ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.snp.vcf --compare ${BASE}/samples/${SAMPLES}/${COMPARE} --column ${COLUMN} --fasta ${BASE_REFERENCE}/bwa/supercontigs_combined_formatted_gatk.fa --gtf ${BASE_REFERENCE}/${REFERENCE}-combined_gatk.gtf > ${BASE}/pipeline-runs/${RUN_NAME}/gatk_features/variants.snp_assoc_${NAME}.txt

