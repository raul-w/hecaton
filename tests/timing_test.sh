#!/bin/bash

set -eu

mkdir -p test_time_wd && cd test_time_wd && \

nextflow run ../nextflow/hecaton_timing.nf -c ../nextflow/nextflow.config -resume \
	--genome_file ../tests/functional/test.fa \
	--reads "../tests/functional/test{1,2}.fq.gz" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_time_o  