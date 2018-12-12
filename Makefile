all: samples reference conifer features

clean: samples_clean reference_clean conifer_clean gatk_clean

samples: FORCE
	./scripts/generate_trimmed_fastq_samples.sh

samples_clean: 
	./scripts/clean_samples.sh

reference: FORCE
	./scripts/generate_reference_conifer.sh
	./scripts/generate_reference_gatk.sh

reference_clean:
	./scripts/clean_reference.sh

conifer: FORCE
	./scripts/run_conifer.sh

conifer_clean:
	./scripts/clean_run.sh

features: FORCE
	./scripts/run_feature_selection.sh ${ARGS}

gatk: FORCE
	./scripts/run_gatk.sh

gatk_features: FORCE
	./scripts/run_gatk_assoc.sh ${ARGS}

gatk_clean:
	./scripts/clean_gatk.sh

lamplink: FORCE
	./scripts/run_lamplink.sh ${ARGS}

FORCE:

