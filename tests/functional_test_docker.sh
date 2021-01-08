#!/bin/bash

set -e

mkdir -p test_wd_docker && cd test_wd_docker && \

nextflow run -w functional_workdir -with-docker hecaton:v0.4.0 -c ../nextflow/nextflow.config -resume ../nextflow/hecaton.nf \
	--genome_file ../tests/functional/test.fa \
	--reads "../tests/functional/test{1,2}.fq.gz" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_o && \
	echo "Functional test with non-empty input successful" && \

nextflow run -w functional_workdir_empty -with-docker hecaton:v0.4.0 -c ../nextflow/nextflow.config -resume ../nextflow/hecaton.nf \
	--genome_file ../tests/functional/test.fa \
	--reads "../tests/functional/test_empty{1,2}.fq.gz" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_o_empty && \
	echo "Functional test with empty input successful" && \

nextflow run -w functional_workdir_no_align -with-docker hecaton:v0.4.0 -c ../nextflow/nextflow.config -resume ../nextflow/hecaton_no_align.nf \
	--genome_file ../tests/functional/test.fa \
	--bwa_bams "../tests/functional/test.bam" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_no_align && \
	echo "Functional test without doing alignments successful"