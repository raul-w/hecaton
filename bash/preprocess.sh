#!/bin/bash

set -e

# Index FASTA file
samtools faidx "$1" &&
# Create BWA index
bwa index "$1" &&
# Create BED file with N regions
seqtk cutN -gp10000000 -n1 "$1" > "$1".N.bed
# create .genome file
cut -f 1-2 "$1".fai > "$1".genome
# create dict file
java -jar $PICARD CreateSequenceDictionary R="$1" O="$1".dict

