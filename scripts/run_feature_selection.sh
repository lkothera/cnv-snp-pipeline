#!/bin/bash

# select features from cnv data based on sample labels


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
module load Python/2.7.12


# create features folder
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/features

# label samples
${BASE}/scripts/extract_label_samples.py --data_path ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_features.txt --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column ${COLUMN} > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_values_${NAME}_labeled.txt

${BASE}/scripts/extract_label_samples.py --data_path ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_gene_features.txt --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column ${COLUMN} > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_values_${NAME}_labeled.txt

# select features
${BASE}/bin/mosquito_feature_selection -data_path ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_values_${NAME}_labeled.txt > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_features_${NAME}.txt

${BASE}/bin/mosquito_feature_selection -data_path ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_values_${NAME}_labeled.txt > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_features_${NAME}.txt
