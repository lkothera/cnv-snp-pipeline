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
module load lamplink/1.12
module load vcftools/0.1.14

# create lamplink folder
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/lamplink

# create plink files
vcftools --vcf ${BASE}/pipeline-runs/${RUN_NAME}/gatk/variants.snp.vcf --out ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/plink_coords --plink

# remove extra columns from plink files
cat ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/plink_coords.ped | cut -f1,2,7- > ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/plink_coords.ped2
mv ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/plink_coords.ped2 ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/plink_coords.ped

# create pheno file
${BASE}/scripts/create_lamplink_pheno.py --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column ${COLUMN} > ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/${NAME}_pheno.txt

# run lamplink
lamplink --file ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/plink_coords --lamp --model-rec --sglev 0.2 --upper 0.5 --fisher --no-parents --no-sex --no-pheno --pheno ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/${NAME}_pheno.txt --out ${BASE}/pipeline-runs/${RUN_NAME}/lamplink/${NAME}


