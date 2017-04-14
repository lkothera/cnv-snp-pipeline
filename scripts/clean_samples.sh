#!/bin/bash

# clean up generated sample QC files

# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# exit if variables weren't set properly
if [[ -z "${FAQCS_DIR}" ]]; then
	exit 1
fi

rm -rf ${BASE}/samples/${SAMPLES}/${FAQCS_DIR}
rm -rf ${BASE}/samples/${SAMPLES}/${FAQCS_DIR}_only

rm ${BASE}/samples/${SAMPLES}/sample_names.txt

