#!/usr/bin/env nextflow
/*
 *Run the filtering on size and quality score step of Hecaton
 *
 *@authors
 *Raúl Wijfjes <raul.wijfjes@wur.nl>
 *
 *Date last modified: 29-11-2019
 */

params.probability_files = ""
params.cutoff = 0.7
params.output_dir = ""

/*
 * Input parameters validation and setup
 */

if (! params.output_dir ) {
    println "Missing output directory parameter"
    exit 1
}
if (! params.probability_files ) {
    println "Missing probability files parameter"
    exit 1
}
probability_files = Channel.fromPath(params.probability_files)
			     .map { file -> tuple(getSampleId(file), file)}
			     .groupTuple()

def getSampleId(file) {
	file.getBaseName()
}

/*
 * Filter calls according to size and cutoff
 */

process filter_calls_cutoff {
	publishDir "${params.output_dir}", mode: 'move' 
	tag "Filter calls of ${prefix}. Cutoff: ${params.cutoff}" 

	input:
	set val(prefix), file(probability_file) from probability_files

	output:
	set val(prefix), file("${prefix}_filtered_cutoff_${params.cutoff}.bedpe")

	script:
	"""
	source activate hecaton_py3
    filter_calls_by_query.py \
    -i ${probability_file} \
    -q \"SIZE >= 50 & PREDICTION_1 >= ${params.cutoff}\" \
    -o ${prefix}_filtered_cutoff_${params.cutoff}.bedpe
    source deactivate
	"""
}
