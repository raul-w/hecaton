#!/bin/bash

set -eu

mkdir -p test_wd && cd test_wd && \

merge_vcf_files.py -i ../tests/genotyping/test_input.txt -f ../tests/genotyping/test.fa.fai -o test_output.vcf

DIFF=$(diff test_output.vcf ../tests/genotyping/expected_output.vcf) 
if [ "$DIFF" != "" ] 
then
    echo "Test failed"
else
	echo "Test was successful"
fi  