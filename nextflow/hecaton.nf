#!/usr/bin/env nextflow
/*
 *Run Hecaton
 *
 *@authors
 *Raúl Wijfjes <raul.wijfjes@wur.nl>
 *
 */

params.genome_file = ""
params.reads = "prefix_{1,2}.fastq"
params.alignment_exclude = 'NO_EXCLUDE'
params.model_file = ""
params.cutoff = 0.7
params.extra_filtering = false
params.manta_config = ""
params.output_dir = ""
params.help = false

/*
 * Add help message
 */

def helpMessage() {
    log.info"""
    =========================================
     Hecaton v0.5.2
    =========================================
    Usage:
    nextflow run hecaton --genome_file reference.fa --reads "prefix_{1,2}.fastq" --manta_config configManta_weight_1.py.ini --model_file model_file.pkl --output_dir results

    Mandatory arguments:
    --genome_file: Reference genome (processed by preprocess.sh) in FASTA format
    --reads: Glob pattern specifying the location of a set of paired-end reads in FASTQ format. Gzipped FASTQ files are allowed
    --manta_config: Config file that will be passed to the Manta tool. Can be found in the "docker" directory of the Hecaton repository
    --model_files: Random forest model that will be used to filter reads. Models can be found in the "models" directory of the Hecaton repository
    --output_dir: Output directory to which all results will be written
    
    Optional arguments:
    --alignment_exclude: BED file containing chromosomal regions that will be excluded from CNV calling
    --extra_filtering: If true, filter CNVs based on read depth and overlap with Ns. Default: false 
    -w: Working directory to which intermediate results will be written. Default: work
    -c: Config file specifying the number of CPU cores and memory that will be assigned to Hecaton
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage

if (params.help){
    helpMessage()
    exit 0
}

/*
 * Input parameters validation and setup
 */

if (! params.output_dir ) {
    println "Missing output directory parameter"
    exit 1
}
if (! params.reads ) {
    println "Missing reads parameter"
    exit 1
}
if (! params.genome_file ) {
    println "Missing genome file parameter"
    exit 1
}
read_files = Channel.fromFilePairs(params.reads)
genome_file = file(params.genome_file)
genome_index_file = file(params.genome_file + ".fai")
genome_bwa_amb_file = file(params.genome_file + ".amb")
genome_bwa_ann_file = file(params.genome_file + ".ann")
genome_bwa_bwt_file = file(params.genome_file + ".bwt")
genome_bwa_pac_file = file(params.genome_file + ".pac")
genome_bwa_sa_file = file(params.genome_file + ".sa")
genome_N_file = file(params.genome_file + ".N.bed")
genome_bed = file(params.genome_file + ".genome")
manta_config_file = file(params.manta_config)
model_file = file(params.model_file)
exclude_file = file(params.alignment_exclude)

/*
 * Align reads to reference genome using speedseq
 */

process align_reads {
	label "bwa" 
	publishDir "${params.output_dir}/aligned_reads", mode: 'copy'
	tag "Aligning reads with bwa mem: ${reads[0]} ${reads[1]}"

	input:
	set val(prefix), file(reads) from read_files
	file genome_file
	file genome_bwa_amb_file
	file genome_bwa_ann_file
	file genome_bwa_bwt_file
	file genome_bwa_pac_file
	file genome_bwa_sa_file
	file exclude_file

	output:
	set val(prefix), file("${prefix}.bam") into bwa_bams
	set val(prefix), file("${prefix}.bam.bai") into bwa_bam_indices 
	set val(prefix), file("${prefix}.discordants.bam") into bwa_discordants_bams
	set val(prefix), file("${prefix}.discordants.bam.bai") into bwa_discordants_bam_indices
	set val(prefix), file("${prefix}.splitters.bam") into bwa_splitters_bams
	set val(prefix), file("${prefix}.splitters.bam.bai") into bwa_splitters_bam_indices

	script:
	if( exclude_file.name == 'NO_EXCLUDE' ) 
		"""
		conda activate hecaton_py3 && \
		speedseq align -t ${task.cpus} -o ${prefix} -R \"@RG\tID:${prefix}\tSM:${prefix}\tLB:${prefix}\" \
		${genome_file} ${reads[0]} ${reads[1]} && \
		conda deactivate
		"""
	else
		"""
		conda activate hecaton_py3 && \
		speedseq align -t ${task.cpus} -o ${prefix} -R \"@RG\tID:${prefix}\tSM:${prefix}\tLB:${prefix}\" \
		${genome_file} ${reads[0]} ${reads[1]} && \
		samtools view -h -b -L ${exclude_file} ${prefix}.bam > ${prefix}_excluded.bam && \
		samtools index ${prefix}_excluded.bam && \
		mv ${prefix}_excluded.bam ${prefix}.bam && \
		mv ${prefix}_excluded.bam.bai ${prefix}.bam.bai && \
		samtools view -h -b -L ${exclude_file} ${prefix}.discordants.bam > ${prefix}.discordants_excluded.bam && \
		samtools index ${prefix}.discordants_excluded.bam && \
		mv ${prefix}.discordants_excluded.bam ${prefix}.discordants.bam && \
		mv ${prefix}.discordants_excluded.bam.bai ${prefix}.discordants.bam.bai && \
		samtools view -h -b -L ${exclude_file} ${prefix}.splitters.bam > ${prefix}.splitters_excluded.bam && \
		samtools index ${prefix}.splitters_excluded.bam && \
		mv ${prefix}.splitters_excluded.bam ${prefix}.splitters.bam && \
		mv ${prefix}.splitters_excluded.bam.bai ${prefix}.splitters.bam.bai && \
		conda deactivate
		"""
}

/*
 * Create channels for all callers
 */

bwa_bams.into{delly_bams; gridss_bams; lumpy_bams; manta_bams; duphold_bams}
bwa_bam_indices.into{delly_bam_indices; gridss_bam_indices; lumpy_bam_indices; manta_bam_indices; duphold_bam_indices}


/*
 * Call CNVs with LUMPY
 */

process call_lumpy {
	label "multithreading"
	publishDir "${params.output_dir}/lumpy_calls", mode: 'copy'
	tag "LUMPY calling: ${prefix}"

	input:
	set val(prefix), file(alignment_file) from lumpy_bams
	set val(prefix_index), file(alignment_file_index) from lumpy_bam_indices
	set val(discordants_prefix), file(discordants_alignment_file) from bwa_discordants_bams
	set val(discordants_prefix_index), file(discordants_alignment_file_index) from bwa_discordants_bam_indices
	set val(splitters_prefix), file(splitters_alignment_file) from bwa_splitters_bams
	set val(splitters_prefix_index), file(splitters_alignment_file_index) from bwa_splitters_bam_indices
	file genome_file
	file genome_index_file
	file genome_N_file

	output:
	file "${prefix}.sv.vcf.gz"
	file "${prefix}.sv.vcf.gz.tbi"
	file "${prefix}.bedpe"
	file "${prefix}_collapsed.bedpe"
	file "${prefix}_post_processed.bedpe"
	set val(prefix), file ("${prefix}_post_processed_collapsed.bedpe") into lumpy_calls

	script:
	"""
	conda activate hecaton_py3 &&
	speedseq sv -t ${task.cpus} -o ${prefix} \
	-x ${genome_N_file} \
	-B ${alignment_file} \
	-D ${discordants_alignment_file} \
	-S ${splitters_alignment_file} \
	-R ${genome_file} \
	-m 1 &&
	vcf_to_bedpe.py -i ${prefix}.sv.vcf.gz \
	-o ${prefix}.bedpe -t LUMPY &&
	process_simple_cnvs.py -i ${prefix}.sv.vcf.gz -o ${prefix}_post_processed.bedpe -t LUMPY &&
	internally_collapse_bedpe_intervals.py \
	-b ${prefix}.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_collapsed.bedpe
	internally_collapse_bedpe_intervals.py -e True \
	-b ${prefix}_post_processed.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_post_processed_collapsed.bedpe
	conda deactivate
	"""
}

/*
 * Call CNVs with DELLY
 */

process call_delly {
	publishDir "${params.output_dir}/delly_calls", mode: 'copy'
	tag "DELLY calling: ${prefix}"

	input:
	set val(prefix), file(alignment_file) from delly_bams
	set val(prefix_index), file(alignment_file_index) from delly_bam_indices
	file genome_file
	file genome_index_file

	output:
	file "${prefix}_delly_cnvs.vcf.gz"
	file "${prefix}_delly_cnvs.vcf.gz.tbi"
	file "${prefix}.bedpe"
	file "${prefix}_collapsed.bedpe"
	file "${prefix}_post_processed.bedpe"
	set val(prefix), file("${prefix}_post_processed_collapsed.bedpe") into delly_calls

	script:
	"""
	conda activate hecaton_py3
	delly call -g ${genome_file} -o ${prefix}_delly_cnvs.bcf ${alignment_file} &&
	bcftools view ${prefix}_delly_cnvs.bcf > ${prefix}_delly_cnvs.vcf &&
	bgzip -f ${prefix}_delly_cnvs.vcf &&
	tabix -p vcf ${prefix}_delly_cnvs.vcf.gz
	vcf_to_bedpe.py -i ${prefix}_delly_cnvs.vcf.gz \
	-o ${prefix}.bedpe -t DELLY &&
	process_simple_cnvs.py -i ${prefix}_delly_cnvs.vcf.gz \
	-o ${prefix}_post_processed.bedpe -t DELLY
	internally_collapse_bedpe_intervals.py \
	-b ${prefix}.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_collapsed.bedpe
	internally_collapse_bedpe_intervals.py -e True \
	-b ${prefix}_post_processed.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_post_processed_collapsed.bedpe
	conda deactivate
	"""
}

/*
 * Call CNVs with GRIDSS
 */

process call_gridss {
	label "multithreading"
	publishDir "${params.output_dir}/gridss_calls", mode: 'copy'
	tag "GRIDSS calling: ${prefix}"

	input:
	set val(prefix), file(alignment_file) from gridss_bams
	set val(prefix_index), file(alignment_file_index) from gridss_bam_indices
	file genome_file
	file genome_index_file
	file genome_bwa_amb_file
	file genome_bwa_ann_file
	file genome_bwa_bwt_file
	file genome_bwa_pac_file
	file genome_bwa_sa_file

	output:
	file "${prefix}.sv.vcf"
	file "${prefix}_post_processed.bedpe"
	set val(prefix), file("${prefix}_post_processed_collapsed.bedpe") into gridss_calls

	script:
	"""
	conda activate hecaton_py3
	java -ea -Xmx31g \
	-Dreference_fasta=${genome_file} \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dgridss.gridss.output_to_temp_file=true \
	-Dgridss.keepTempFiles=true \
	-cp \$GRIDSS_JAR gridss.CallVariants \
	TMP_DIR=. \
	WORKING_DIR=. \
	WORKER_THREADS=${task.cpus} \
	REFERENCE_SEQUENCE=${genome_file} \
	INPUT=${alignment_file} \
	OUTPUT=${prefix}.sv.vcf \
	ASSEMBLY=${prefix}.gridss.assembly.bam \
	2>&1 | tee -a gridss.\$HOSTNAME.\$\$.log &&
	fix_gridss_header.py -i ${prefix}.sv.vcf -o gridss_fixed.vcf
	process_simple_cnvs.py -i gridss_fixed.vcf \
	-o ${prefix}_post_processed.bedpe -t GRIDSS &&
	internally_collapse_bedpe_intervals.py -e True \
	-b ${prefix}_post_processed.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_post_processed_collapsed.bedpe
	conda deactivate
	"""
}

/*
 * Call CNVs with Manta
 */

process call_manta {
	label "multithreading"
	publishDir "${params.output_dir}/manta_calls", mode: 'copy'
	tag "Manta calling: ${prefix}"

	input:
	set val(prefix), file(alignment_file) from manta_bams
	set val(prefix_index), file(alignment_file_index) from manta_bam_indices
	file genome_file
	file genome_index_file

	output:
	file "${prefix}_tumorSV.vcf.gz"
	file "${prefix}_tumorSV.vcf.gz.tbi"
	file "${prefix}.bedpe"
	file "${prefix}_collapsed.bedpe"
	file "${prefix}_post_processed.bedpe"
	set val(prefix), file("${prefix}_post_processed_collapsed.bedpe") into manta_calls

	script:
	"""
	conda activate hecaton_py2
	configManta.py --tumorBam ${alignment_file} \
	--config ${manta_config_file} \
	--referenceFasta ${genome_file} --runDir ${prefix}_rundir &&
	cd ${prefix}_rundir &&
	python2.7 runWorkflow.py -m local -j ${task.cpus} &&
	mv results/variants/tumorSV.vcf.gz ../${prefix}_tumorSV.vcf.gz &&
	mv results/variants/tumorSV.vcf.gz.tbi ../${prefix}_tumorSV.vcf.gz.tbi &&
	cd .. &&
	conda deactivate
	conda activate hecaton_py3
	vcf_to_bedpe.py -i ${prefix}_tumorSV.vcf.gz \
	-o ${prefix}.bedpe -t Manta &&
	process_simple_cnvs.py -i ${prefix}_tumorSV.vcf.gz \
	-o ${prefix}_post_processed.bedpe -t Manta
	internally_collapse_bedpe_intervals.py \
	-b ${prefix}.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_collapsed.bedpe
	internally_collapse_bedpe_intervals.py -e True \
	-b ${prefix}_post_processed.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_post_processed_collapsed.bedpe
	conda deactivate
	"""
}

joined_calls = delly_calls.join(gridss_calls).join(lumpy_calls).join(manta_calls)

/*
 * Intersect all calls
 */
 
process intersecting_calls {
	publishDir "${params.output_dir}/intersected_calls", mode: 'copy'
	tag "Intersecting calls: ${prefix}"

	input:
	set val(prefix), file("delly_file.bedpe"), file("gridss_file.bedpe"), file("lumpy_file.bedpe"), file("manta_file.bedpe") from joined_calls
	file genome_index_file

	output:
	set val(prefix), file("${prefix}_intersected.bedpe") into intersected_calls

	script:
	"""
	conda activate hecaton_py3
	intersecting_bedpe_intervals.py \
	-a lumpy_file.bedpe \
	-b manta_file.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_1.bedpe && 
	intersecting_bedpe_intervals.py \
	-a ${prefix}_1.bedpe \
	-b delly_file.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_2.bedpe &&
	intersecting_bedpe_intervals.py \
	-a ${prefix}_2.bedpe \
	-b gridss_file.bedpe \
	-f ${genome_index_file} \
	-o ${prefix}_intersected.bedpe
	conda deactivate
	"""
}

/*
 * Prepare features
 */

process prepare_features {
	publishDir "${params.output_dir}/feature_bedpes", mode: 'copy'
	tag "Feature BEDPE preparation: ${prefix}"

	input:
	set val(prefix), file(intersection_file) from intersected_calls

	output:
	set val(prefix), file("${prefix}_insertion_features.bedpe") into insertion_feature_files
	set val(prefix), file("${prefix}_non_insertion_features.bedpe") into non_insertion_feature_files

	script:
	"""
	conda activate hecaton_py3
	intersect_bedpe_to_feature_bedpe.py -i ${intersection_file} \
	-o ${prefix}_features.bedpe
	split_feature_set_non_insertions.py -f ${prefix}_features.bedpe -i ${prefix}_insertion_features.bedpe -n ${prefix}_non_insertion_features.bedpe
	conda deactivate
	"""
}

/*
 * Predict results using random forest
 */

process apply_random_forest {
	publishDir "${params.output_dir}/random_forest_calls", mode: 'copy' 
	tag "Apply random forest: ${prefix}"

	input:
	set val(prefix), file(insertion_file) from insertion_feature_files
	set val(prefix), file(non_insertion_file) from non_insertion_feature_files
	file model_file

	output:
	set val(prefix), file("${prefix}_probabilities_unfiltered.bedpe") into probability_files

	script:
	"""
	conda activate hecaton_py3
	predict_cnvs_model_insertions_probabilities.py \
	-i ${insertion_file} \
    -m ${model_file} \
    -o ${prefix}_insertion_probabilities.bedpe &&
	predict_cnvs_model_insertions_probabilities.py \
	-i ${non_insertion_file} \
    -m ${model_file} \
    -o ${prefix}_non_insertion_probabilities.bedpe &&
	concatenate_insertion_non_insertion_sets.py \
	-i ${prefix}_insertion_probabilities.bedpe \
	-n ${prefix}_non_insertion_probabilities.bedpe \
	-o ${prefix}_probabilities_unfiltered.bedpe
	conda deactivate
	"""
}

/*
 * Filter calls according to size and cutoff
 */

process filter_calls_cutoff {
	publishDir "${params.output_dir}/random_forest_calls", mode: 'copy' 
	tag "Filter calls of ${prefix}. Cutoff: ${params.cutoff}" 

	input:
	set val(prefix), file(probability_file) from probability_files

	output:
	set val(prefix), file("${prefix}_probabilities_filtered_cutoff_${params.cutoff}.bedpe") into score_files

	script:
	"""
	conda activate hecaton_py3
    filter_calls_by_query.py \
    -i ${probability_file} \
    -q \"SIZE >= 50 & PREDICTION_1 >= ${params.cutoff}\" \
    -o ${prefix}_probabilities_filtered_cutoff_${params.cutoff}.bedpe
    conda deactivate
	"""
}

joined_duphold_bams = score_files.join(duphold_bams).join(duphold_bam_indices)

/*
 * Filter calls according to median depth
 */

process filter_calls_median_depth {
	label "multithreading"
	tag "Filter calls of ${prefix} using median depth"

	input:
	set val(prefix), file(probability_file), file(alignment_file), file(alignment_file_index) from joined_duphold_bams
	file genome_file
	file genome_index_file

	output:
	set val(prefix), file("${prefix}_probabilities_filtered_cutoff_${params.cutoff}_depth.bedpe") into depth_files

	when:
	params.extra_filtering != false

	script:
	"""
	conda activate hecaton_py3
	bedpe_to_vcf.py \
	-i ${probability_file} \
	-o ${prefix}.vcf \
	-s ${prefix} &&
	duphold -t ${task.cpus} \
	-v ${prefix}.vcf \
	-b ${alignment_file} \
	-f ${genome_file} \
	-o ${prefix}_duphold.vcf &&
	concat_duphold_to_calls.py \
	-i ${probability_file} \
	-d ${prefix}_duphold.vcf \
	-o ${prefix}_probabilities_duphold.bedpe &&
	filter_calls_by_query.py \
	-i ${prefix}_probabilities_duphold.bedpe \
	-q \"(Chrom_norm_depth < 4 & GC_norm_depth < 4 & Flank_norm_depth < 4) & ((Chrom_norm_depth < 0.7 & GC_norm_depth < 0.7 & Flank_norm_depth < 0.7 & TYPE == 'DEL') | (Chrom_norm_depth > 1.3 & GC_norm_depth > 1.3 & Flank_norm_depth > 1.3 & (TYPE == 'DUP:TANDEM' | TYPE == 'DUP:DISPERSED')) | TYPE == 'INS')\" \
	-o ${prefix}_probabilities_filtered_cutoff_${params.cutoff}_depth.bedpe
	conda deactivate
	"""
}

/*
 * Filter calls according to median depth
 */

process filter_calls_flanking_Ns {
	publishDir "${params.output_dir}/random_forest_calls", mode: 'move'
	tag "Filter calls of ${prefix} using flanking Ns" 

	input:
	set val(prefix), file(probability_file) from depth_files
	file genome_bed

	output:
	set val(prefix), file("${prefix}_probabilities_filtered_cutoff_${params.cutoff}_depth_flanking_Ns.bedpe") into flank_files

	when:
	params.extra_filtering != false

	script:
	"""
	conda activate hecaton_py3
	breakpoint_slop.py -i ${probability_file} \
	-o tmp_${prefix}_slop_200.bedpe -s 200 && \
	bedtools flank -i tmp_${prefix}_slop_200.bedpe \
	-g ${genome_bed} -b 400 \
	> tmp_${prefix}_flanking_400_bp.bed && \
	bedtools getfasta -fi ${genome_file} \
	-bed tmp_${prefix}_flanking_400_bp.bed \
	> tmp_${prefix}_flanking_400_bp.fa && \
	calculate_N_fraction.py -i tmp_${prefix}_flanking_400_bp.fa \
	-o tmp_${prefix}_fraction_Ns_flanking_400_bp.csv && \
	concat_fraction_n_to_calls.py -i ${probability_file} \
	-f tmp_${prefix}_fraction_Ns_flanking_400_bp.csv \
	-o ${prefix}_probabilities_flanking_Ns.bedpe
	filter_calls_by_query.py \
	-i ${prefix}_probabilities_flanking_Ns.bedpe \
	-q \"Flanking_Ns < 0.1\" \
	-o ${prefix}_probabilities_filtered_cutoff_${params.cutoff}_depth_flanking_Ns.bedpe
	conda deactivate
	"""
}

