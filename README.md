# Hecaton: reliable detection of CNV in plant genomes

Hecaton is a framework specifically designed for plant genomes that detects copy number variants (CNVs) using short paired-end Illumina reads. CNVs are called by integrating existing structural variant callers through a machine-learning model and several custom post-processing scripts.

## Table of Contents

- [Workflow of Hecaton](#overview)
- [Installation](#install)
  - [Prerequisites](#prerequisites)
  - [Local installation](#local)
  - [Docker](#docker)
- [Usage](#usage)
  - [General usage](#general-usage)
  - [Running multiple samples at once](#multiple)
  - [Balancing sensitivity and precision](#sensitivity-precision)
  - [Applying Hecaton to low quality reference genomes or samples distantly related to the reference](#complex_genomes)
  - [Converting the BEDPE output to VCF](#convert_to_vcf)
  - [Merging and genotyping multiple samples](#merging)
- [Limitations](#limit)

## <a name="overview"></a>Workflow of Hecaton

Hecaton consists of the following steps:

* Aligning reads to a reference genome using [bwa mem][bwa mem]
* Calling CNVs using the structural variant callers [Delly][delly], [GRIDSS][gridss], [LUMPY][lumpy], and [Manta][manta]
* Post-processing each set of CNVs to remove false positives
* Merging all sets of CNVs into one large set
* Classifying CNVs in this large set as true or false positives using a random forest model
* Optional filtering of CNVs based on read depth (computed by [duphold][duphold]) and the presence of nearby gaps

## <a name="install"></a>Installation

### <a name="prerequisites"></a>Prerequisites

Hecaton can only be installed on Linux systems and requires `git` and Anaconda/Miniconda (Python 3.6+) to be present on the system. Furthermore, ensure libncurses5-dev is installed:
```bash
sudo apt update
sudo apt install libncurses5-dev
```

Additional notes regarding hardware requirements:

* **Memory**: The minimum memory requirements of Hecaton are at least 4 GB of memory per core.
* **CPU**: Only 1 CPU core is required to run Hecaton, but more will generally decrease running time.
* **I/O**: Hecaton may perform a large amount of reading and writing to disc, so increasing disc access speeds by using for instance a local cache can greatly improve running time. 

All steps of Hecaton are run using the Nextflow workflow language, which can be obtained by running: 
```bash
wget -qO- https://get.nextflow.io | bash
```

In order to run Hecaton, the `nextflow` binary needs to be added to the $PATH environment variable (e.g. export PATH=$PATH:directory/to/nextflow).

Finally, the Hecaton scripts themselves need to be obtained and added to $PATH.

Clone the repository to your directory of choice:
```bash
git clone https://git.wageningenur.nl/wijfj001/hecaton.git
```

Set permissions and add all scripts of Hecaton to $PATH:
```bash
cd hecaton
chmod +x scripts/collapse/* && \
chmod +x scripts/convert/* && \
chmod +x scripts/filter/* && \
chmod +x scripts/genotype/* && \  
chmod +x scripts/gridss/* && \
chmod +x scripts/intersect/* && \
chmod +x scripts/predict/* && \
chmod +x scripts/process/* && \
export PATH=$PWD/scripts/collapse:$PWD/scripts/convert:$PWD/scripts/filter:$PWD/scripts/genotype:$PWD/scripts/gridss:$PWD/scripts/intersect:$PWD/scripts/predict:$PWD/scripts/process:$PATH && \
export PYTHONPATH=$PYTHONPATH:$PWD/scripts
```

The dependencies of Hecaton can be installed locally or obtained through a Docker image.

### <a name="local"></a>Local installation

After installing all the prerequisites, the following commands can be used to do a quick local installation in which all dependencies are located in the directory of the Hecaton repository. Note that this script adds a few lines to `~/.bashrc`.

```bash
bash install.sh
``` 

The following command tests whether Hecaton was installed correctly. Please ensure that you have installed Nextflow before running it (see [Prerequisites](#prerequisites)).
```bash
bash tests/functional_test.sh
```
A complete overview of all installation commands is given below.

Most dependencies can be installed by generating two custom `conda` environments, named `hecaton_py3` and `hecaton_py2`:
```bash
cd hecaton && \
conda env create -f docker/environment_py3.yml && \
conda env create -f docker/environment_py2.yml
``` 

The remaining dependencies (gridss, picard, speedseq) need to be installed separately in a directory of choice.

```bash
mkdir hecaton_deps && \
cd hecaton_deps && \
wget https://github.com/PapenfussLab/gridss/releases/download/v2.0.1/gridss-2.0.1-gridss-jar-with-dependencies.jar && \
export GRIDSS_JAR=$PWD/gridss-2.0.1-gridss-jar-with-dependencies.jar && \
wget https://github.com/broadinstitute/picard/releases/download/2.18.23/picard.jar && \
export PICARD=$PWD/picard.jar && \
source activate hecaton_py2 && \
git clone --recursive https://github.com/hall-lab/speedseq && \
cd speedseq && \
make align && \
make sv && \
make config && \
export PATH=$PWD/bin:$PATH && \
source deactivate && \
cd ../.. 
```

The installation directory (hecaton_deps) can be modified according to preference, as long as all binaries in `speedseq/bin` are added to the $PATH environment variable, and the $GRIDSS_JAR and $PICARD environment variables contain the jar files of GRIDSS and Picard respectively. 

### <a name="docker"></a>Docker

Hecaton can be run through a Docker image, avoiding the need to install dependencies locally:
```bash
cd hecaton/docker && \
bash docker_build.sh && \
cd ..
``` 

These commands need to be re-run each time a new version of Hecaton is acquired. 

The following command tests whether the image has been correctly built. Please ensure that you have installed Nextflow before running it (see [Prerequisites](#prerequisites)).

```bash
bash tests/functional_test_docker.sh
```

## <a name="usage"></a>Usage

### <a name="general-usage"></a>General usage

Before running Hecaton, ensure that:
* the `nextflow` binary is added to the $PATH environment variable.
* the $GRIDSS_JAR and $PICARD environment variables contain the jar files of GRIDSS and Picard respectively 
* all binaries in `speedseq/bin` are added to $PATH
* all Hecaton scripts are added to $PATH
* the directory containing the Hecaton scripts is added to $PYTHONPATH

See [Installation](#install) for examples of commands that can be executed to achieve all of the above.

Hecaton takes a reference genome and a set of paired-end reads as input, producing a set of CNVs in [BEDPE][bedpe] format as output. The reference genome first needs to be preprocessed:
```bash
bash bash/preprocess.sh genome.fa
```

This step only needs to be run once for every reference genome. 

A single command can then be used to run Hecaton on a set of paired-end reads:

Local:
```bash
nextflow run -c nextflow/nextflow.config -w hecaton_workdir nextflow/hecaton.nf --genome_file genome.fa --reads "reads_{1,2}.fq" --manta_config docker/configManta_weight_1.py.ini --output_dir output --model models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl
```

Docker:
```bash
nextflow run -with-docker hecaton:v1 -c nextflow/nextflow.config -w hecaton_workdir nextflow/hecaton.nf --genome_file genome.fa --reads "reads_{1,2}.fq" --manta_config docker/configManta_weight_1.py.ini --output_dir output --model models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl
```

The following parameters need to be specified when running Hecaton:
* `--genome_file`: Reference genome (processed by preprocess.sh) in FASTA format
* `--reads`: Glob pattern specifying the location of a set of paired-end reads in FASTQ format. Hecaton can work with gzipped FASTQ files as input, avoiding the need to decompress read data. 
* `--manta_config`: Config file that will be passed to the Manta tool. Can be found in `docker`.
* `--output_dir`: Output directory to which all results will be written
* `--model_file`: Random forest model that will be used to filter reads. Models can be found in the `models` directory of the Hecaton repository. 

The other parameters are optional:
* `-w`:the working directory to which intermediate results will be written (`work` by default). If an instance of Hecaton fails mid-run, it can be resumed from the point of failure by adding the `-resume` parameter to the command, assuming that the command is run in the directory in which the working directory is present:
```bash
nextflow run -c nextflow/nextflow.config -w hecaton_workdir nextflow/hecaton.nf -resume --genome_file genome.fa --reads "reads_{1,2}.fq" --manta_config docker/configManta_weight_1.py.ini --output_dir output --model models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl
```

* `-c`: indicates a config file specifying the number of CPU cores and memory that will be assigned to Hecaton. The default `nextflow/nextflow.config` file assigns a maximum of 16 cores and 32 GB of RAM to Hecaton. This config file can be tuned to fit a specific computational setup. See the [Nextflow manual][nextflow_man] for additional details about generating a config file. As an example, the `nextflow/nextflow_slurm.config` can be used to run Hecaton on a HPC with Slurm (queue stills needs to be specified in this file).    

### <a name="multiple"></a>Running multiple samples at once

A glob pattern can be used to run multiple samples using a single command. If there for instance two sets of paired-end reads named sample1_1.fq, sample1_2.fq, sample2_1.fq, sample2_2.fq, they can all be processed using:
```bash
nextflow run -c nextflow/nextflow.config -w hecaton_workdir nextflow/hecaton.nf --genome_file genome.fa --reads "*_{1,2}.fastq" --manta_config docker/configManta_weight_1.py.ini --output_dir output --model models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl
```

### <a name="sensitivity-precision"></a>Balancing sensitivity and precision

The default cutoff used by the random forest model of Hecaton (0.7) resulted in a good balance of sensitivity and precision in our benchmarks: Hecaton attained at least 80 % precision for all the different types of CNVs. The cutoff can be changed through the `--cutoff` parameter to make Hecaton more stringent or lenient:
```bash
nextflow run -c nextflow/nextflow.config -w hecaton_workdir nextflow/hecaton.nf --genome_file genome.fa --reads "reads_{1,2}.fastq" --manta_config docker/configManta_weight_1.py.ini --output_dir output --cutoff 0.5 --model models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl
```

### <a name="complex_genomes"></a>Applying Hecaton to low quality reference genomes or samples distant from the reference

Low quality reference genomes containing a large number of gaps or samples distantly related to the reference may introduce additional false positives that are difficult to catch using the random forest model. To increase precision for such use cases, Hecaton can be run with optional filtering steps (turned off by default) that remove putative false positive CNVs based on read depth and the presence of nearby stretches of N's in the reference genome. These steps can be invoked by setting the `--extra_filtering` parameter to `true`:
```bash
nextflow run -c nextflow/nextflow.config -w hecaton_workdir nextflow/hecaton.nf --genome_file genome.fa --reads "reads_{1,2}.fastq" --manta_config docker/configManta_weight_1.py.ini --output_dir output --model models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl --extra_filtering true
```

### <a name="convert_to_vcf"></a>Converting BEDPE output to VCF

The `bedpe_to_vcf.py` script can be used to convert BEDPE output files to VCF format. It is recommended to compress and index the VCF file using `bgzip` and `tabix`, and to compute the read depth of each CNV using `duphold`.
```bash
source activate hecaton_py3
scripts/convert/bedpe_to_vcf.py -i output.bedpe -o output.vcf -s name_of_your_sample
bgzip output.vcf
tabix output.vcf.gz
duphold -t number_of_threads -v output.vcf.gz -b alignment_of_this_sample.bam -f reference.fa -o output_duphold.vcf
bgzip output_duphold.vcf
tabix output_duphold.vcf.gz
source deactivate
```

### <a name="merging"></a>Merging and genotyping multiple samples

The output of multiple samples can be merged into a large VCF file, if they have been converted from BEDPE to VCF format and post-processed by `duphold` (see [Converting BEDPE output to VCF](#convert_to_vcf)). Calls are collapsed into a single events if share a specified fraction of reciprocal overlap (default 0.5) and if their breakpoints are located within 1000 bp of each other. 
```bash
source activate hecaton_py3
scripts/genotype/merge_vcf_files.py -i samples.txt -f reference.fa.fai -o merged_samples.vcf -r 0.5
source deactivate
``` 

The merging script uses the following arguments:
* `-i`: File containing a path to a single sample VCF file (compressed by `bgzip`, indexed by `tabix`, and post-processed by `duphold`) on each line
* `-f`: Index file of used reference genome (generated using `samtools faidx`)
* `-o`: Name of the output file.
* `-r`: Minimum fraction of reciprocal overlap at which calls are merged

## <a name="limit"></a>Limitations

* Hecaton can only work with Illumina paired-end reads. Forward and reverse reads need to be present in separate FASTQ files. All other types of sequencing format are not supported.

* Hecaton is designed to detect CNVs: deletions, insertions, tandem duplications, and dispersed duplications (both inter- and intrachromosomal). It does not report other types of structural variants. 

* We were unable to evaluate the performance of Hecaton with regards to detecting CNVs larger than 1 Mb. We recommended excluding such CNVs from further downstream analysis. 

[bwa mem]: https://github.com/lh3/bwa
[delly]: https://github.com/dellytools/delly
[gridss]: https://github.com/PapenfussLab/GRIDSS
[lumpy]: https://github.com/arq5x/lumpy-sv
[manta]: https://github.com/Illumina/manta
[duphold]: https://github.com/brentp/duphold
[bedpe]: https://bedtools.readthedocs.io/en/latest/content/general-usage.html
[nextflow_man]: https://www.nextflow.io/docs/latest/index.html