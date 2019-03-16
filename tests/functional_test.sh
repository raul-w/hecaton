#!/bin/bash

set -eu

mkdir -p test_wd && cd test_wd && \

nextflow run -w functional_workdir -c ../nextflow/nextflow.config -resume ../nextflow/hecaton.nf \
	--genome_file ../tests/functional/test.fa \
	--reads "../tests/functional/test{1,2}.fq.gz" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_o && \
	echo "Functional test with non-empty input successful" && \

nextflow run -w functional_workdir_empty -c ../nextflow/nextflow.config -resume ../nextflow/hecaton.nf \
	--genome_file ../tests/functional/test.fa \
	--reads "../tests/functional/test_empty{1,2}.fq.gz" \
	--model_file ../models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
	--cutoff 0.7 \
	--manta_config ../docker/configManta_weight_1.py.ini \
	--extra_filtering true \
	--output_dir test_o_empty && \
	echo "Functional test with empty input successful" && \

merge_vcf_files.py -i ../tests/genotyping/test_input.txt -f ../tests/genotyping/test.fa.fai -o test_output.vcf

DIFF=$(diff test_output.vcf ../tests/genotyping/expected_output.vcf | wc -l)
if [ "$DIFF" != "4" ] 
	then
	echo "Collapsing VCFs test failed"
    exit 1
else
	echo "Collapsing VCFs test successful"
fi

vcf_to_table.py -v test_output.vcf -o test_output.tsv -f CHROM POS REF ALT HOM-VAR VAR -gf RQ

cmp test_output.tsv ../tests/genotyping/expected_output.tsv
if [ $? -ne 0 ] 
	then
    echo "Converting VCF to table test failed"
    exit 1
else
	echo "Converting VCF to table test successful"
fi

filter_ref_sites_vcf.py -v ../tests/genotyping/test_filter_ref_input.vcf -o test_filter_ref_output.vcf

cmp test_filter_ref_output.vcf ../tests/genotyping/test_filter_ref_expected_output.vcf
if [ $? -ne 0 ] 
	then
    echo "Filtering reference variants from VCF test failed"
    exit 1
else
	echo "Filtering reference variants from VCF test successful"
fi