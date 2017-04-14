#!/bin/bash

# select features from cnv data based on sample labels


# load pipeline variables
source $(dirname "$0")/pipeline.sh

# load dependencies
source /etc/profile.d/modules.sh
module load Python/2.7.12


# create features folder
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/features

# label samples
${BASE}/scripts/extract_label_samples.py --data_path ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_features.txt --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column 1 > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_values_malathion_labeled.txt

${BASE}/scripts/extract_label_samples.py --data_path ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_features.txt --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column 2 > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_values_permethrin_labeled.txt

${BASE}/scripts/extract_label_samples.py --data_path ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_gene_features.txt --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column 1 > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_values_malathion_labeled.txt

${BASE}/scripts/extract_label_samples.py --data_path ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_gene_features.txt --label_path ${BASE}/samples/${SAMPLES}/${COMPARE} --column 2 > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_values_permethrin_labeled.txt

# select features
${BASE}/bin/mosquito_feature_selection -data_path ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_values_malathion_labeled.txt > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_features_malathion.txt

${BASE}/bin/mosquito_feature_selection -data_path ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_values_permethrin_labeled.txt > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_features_permethrin.txt

${BASE}/bin/mosquito_feature_selection -data_path ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_values_malathion_labeled.txt > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_features_malathion.txt

${BASE}/bin/mosquito_feature_selection -data_path ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_values_permethrin_labeled.txt > ${BASE}/pipeline-runs/${RUN_NAME}/features/cnv_gene_features_permethrin.txt
