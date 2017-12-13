#!/bin/bash

# calculate cnv using conifer


# load pipeline variables
dos2unix $(dirname "$0")/pipeline.sh
source $(dirname "$0")/pipeline.sh

# load dependencies
source /etc/profile.d/modules.sh
module load Python/2.7.12
module load bwa/0.7.12
module load samtools/1.3.1
module load BEDTools/2.17.0
module load R/3.2.3
module load conifer/0.2.2-sing

# create run directory
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}

# align samples to reference
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/bam
${BASE}/scripts/bwa_align.py --in_dir ${BASE}/samples/${SAMPLES}/${FAQCS_DIR}_only --out_dir ${BASE}/pipeline-runs/${RUN_NAME}/bam --ext fastq --database ${BWA_REFERENCE_CONIFER}

# index bam files
ls ${BASE}/pipeline-runs/${RUN_NAME}/bam | grep '.bam$' | sed 's/\(.*\)\.bam/\1/' | xargs --verbose -n 1 -I % sh -c "samtools index ${BASE}/pipeline-runs/${RUN_NAME}/bam/%.bam"

# generate bedgraph file
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/wgx
ls ${BASE}/pipeline-runs/${RUN_NAME}/bam | grep '.bam$' | sed 's/\(.*\)\.bam/\1/' | xargs --verbose -n 1 -I % sh -c "bedtools genomecov -bg -ibam ${BASE}/pipeline-runs/${RUN_NAME}/bam/%.bam > ${BASE}/pipeline-runs/${RUN_NAME}/wgx/%.wgx"

# plot coverage
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/coverage
R --no-save --args ${BASE}/reference/${REFERENCE}/gene_coords_conifer.txt ${BASE}/samples/${SAMPLES}/sample_names.txt ${BASE}/pipeline-runs/${RUN_NAME}/wgx/ ${BASE}/pipeline-runs/${RUN_NAME}/coverage/ < ${BASE}/scripts/generate_coverage.R

## conifer analysis ----
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/conifer

# -- calculate rpkm
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/conifer/rpkm
ls ${BASE}/pipeline-runs/${RUN_NAME}/bam | grep '.bam$' | sed 's/\(.*\)\.bam/\1/' | xargs --verbose -n 1 -I % sh -c "conifer rpkm --probes ${BASE}/reference/${REFERENCE}/probes_final_conifer.txt --input ${BASE}/pipeline-runs/${RUN_NAME}/bam/%.bam --output ${BASE}/pipeline-runs/${RUN_NAME}/conifer/rpkm/%.rpkm.txt"

# -- analyze probes
conifer analyze --probes ${BASE}/reference/${REFERENCE}/probes_final_conifer.txt --rpkm_dir ${BASE}/pipeline-runs/${RUN_NAME}/conifer/rpkm/ --output ${BASE}/pipeline-runs/${RUN_NAME}/conifer/analysis.hdf5 --svd 6 --write_svals ${BASE}/pipeline-runs/${RUN_NAME}/conifer/singular_values.txt --plot_scree ${BASE}/pipeline-runs/${RUN_NAME}/conifer/screeplot.png

# -- cnv calls
conifer call --threshold 0.5 --input ${BASE}/pipeline-runs/${RUN_NAME}/conifer/analysis.hdf5 --output ${BASE}/pipeline-runs/${RUN_NAME}/conifer/cnv_calls.txt

# -- plot calls
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/conifer/call_imgs
conifer plotcalls --input ${BASE}/pipeline-runs/${RUN_NAME}/conifer/analysis.hdf5 --calls ${BASE}/pipeline-runs/${RUN_NAME}/conifer/cnv_calls.txt --outputdir ${BASE}/pipeline-runs/${RUN_NAME}/conifer/call_imgs/

# -- export calls
mkdir -p ${BASE}/pipeline-runs/${RUN_NAME}/conifer/export
conifer export --input ${BASE}/pipeline-runs/${RUN_NAME}/conifer/analysis.hdf5 --output ${BASE}/pipeline-runs/${RUN_NAME}/conifer/export/

# -- combine exports
${BASE}/scripts/combine_exports.py --in_dir ${BASE}/pipeline-runs/${RUN_NAME}/conifer/export --out ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_features.txt

# -- combine exports, combine all exons for each gene
${BASE}/scripts/combine_gene_exports.py --in_dir ${BASE}/pipeline-runs/${RUN_NAME}/conifer/export --out ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_gene_features.txt

# -- combine rpkm values
${BASE}/scripts/combine_rpkms.py --in_dir ${BASE}/pipeline-runs/${RUN_NAME}/conifer/rpkm --combined_exports_file ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_features.txt --probes_file ${BASE}/reference/${REFERENCE}/probes_final_conifer.txt --out ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_rpkm_features.txt --out_gene ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_rpkm_gene_features.txt --log=false

# -- combine rpkm log2 values
${BASE}/scripts/combine_rpkms.py --in_dir ${BASE}/pipeline-runs/${RUN_NAME}/conifer/rpkm --combined_exports_file ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_cnv_features.txt --probes_file ${BASE}/reference/${REFERENCE}/probes_final_conifer.txt --out ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_rpkm_log_features.txt --out_gene ${BASE}/pipeline-runs/${RUN_NAME}/conifer/combined_rpkm_log_gene_features.txt --log=true


