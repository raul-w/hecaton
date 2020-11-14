/*
 *Convert CNV calls to VCF format and add read depth information
 *
 *@authors
 *Raúl Wijfjes <raul.wijfjes@wur.nl>
 *
 *Date last modified: 20-04-2020
 */

def helpMessage() {
    log.info"""

    Heavily inspired by https://github.com/nf-core/methylseq

    Usage:
    An typical command for running this script is as follows:
    nextflow run add_read_depth.nf --fasta genome.fa --input input.tsv -profile conda
    Mandatory arguments:
      --fasta [file]                    Path to fasta reference
      --input [file]                    Path to tab-separated input file of the form "sample\tbam_file\tbedpe_file"
      -profile [str]                    Configuration profile to use. Can use multiple (comma separated)
                                            Available: conda,annuna

    Other options:
     --outdir [file]                    The output directory where the results will be saved
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "fasta file not found : ${params.fasta}" }
        .into { ch_fasta }

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

/*
 * Create a channel for input files
 */


Channel
		.fromPath(params.input, checkIfExists: true)
		.ifEmpty { exit 1, "input file not found : ${params.input}" }
		.splitCsv(sep: '\t')
		.map{ row -> 
				def sampleId = row[0]
				def bamFile = returnFile(row[1])
				def bamIndex = returnFile(row[1] + ".bai")
				def bedpeFile = returnFile(row[2])
				[sampleId, bamFile, bamIndex, bedpeFile] }
		.into { ch_inputFiles }

// Header log info
def summary = [:]
summary['Run Name']  = custom_runName ?: workflow.runName
summary['Reference'] = params.fasta
summary['Input TSV'] = params.input
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Pipeline dir']     = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
if(params.project) summary['Cluster Project'] = params.project
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(params.email)            summary['E-mail Address'] = params.email
if(params.email_on_fail)    summary['E-mail on failure'] = params.email_on_fail
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    file "v_duphold.txt"
    file "v_bgzip.txt"
    file "v_tabix.txt"

    script:
    """
    duphold -h &> v_duphold.txt || true
    bgzip --version &> v_bgzip.txt || true
    tabix --version &> v_tabix.txt || true
    """
}

/*
 * STEP 1 - Convert BEDPEs to VCFs
 */

process convert_bedpe_to_vcf {
	tag "$name"

	input:
	set val(name), file(bam), file(bam_index), file(cnvs) from ch_inputFiles

	output:
	set val(name), file(bam), file(bam_index), file("${name}.vcf.gz"), file("${name}.vcf.gz.tbi") into ch_files_for_duphold

	script:
	"""
	(head -n 1 $cnvs && \\
	tail -n +2 $cnvs | sort -k1,1 -k2,2n) > sorted.bedpe && \\
	bedpe_to_vcf.py \\
	-i sorted.bedpe \\
	-o ${name}.vcf \\
	-s ${name} && \\
	bgzip -f ${name}.vcf && \\
	tabix -p vcf ${name}.vcf.gz
	"""
}

/*
 * STEP 2 - Run duphold on VCFs
 */

 process duphold {
	publishDir "${params.outdir}/duphold", mode: 'copy'
	tag "$name"

	input:
	set val(name), file(bam), file(bam_index), file(vcf), file(vcf_index) from ch_files_for_duphold
	file genome from ch_fasta.collect()

	output:
	file "${name}_duphold.vcf.gz"
	file "${name}_duphold.vcf.gz.tbi"

	script:
	cores = task.cpus
	"""
	mv $bam ${name}.bam
	mv $bam_index ${name}.bam.bai
	duphold -t $cores \\
    -v $vcf \
    -b ${name}.bam \
    -f $genome \
    -o ${name}_duphold.vcf && \\
    bgzip -f ${name}_duphold.vcf && \\
	tabix -p vcf ${name}_duphold.vcf.gz
	"""
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}