#!/bin/bash

# clean up conifer pipeline run


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# exit if variables weren't set properly
if [[ -z "${RUN_NAME}" ]]; then
	exit 1
fi

rm -rf ${BASE}/pipeline-runs/${RUN_NAME}


