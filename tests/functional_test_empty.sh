#!/bin/bash

set -eu

mkdir -p test_wd_empty && cd test_wd_empty && \

nextflow run ../nextflow/hecaton.nf -c ../nextflow/nextflow.config -resume \
	--genome_file ../tests/functional/test.fa \
	--reads "../tests/functional/test_empty{1,2}.fq.gz" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_o  