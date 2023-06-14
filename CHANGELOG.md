# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.2] - 2023-06-14
### Changed
- Added output section to README
- Specified that Nextflow version 22.10.0 should be used, as newer versions do not support DSL 1
- Added "nextflow.enable.dsl=1" to config files, so that DSL 1 will be used instead of the default DSL 2
- Replaced "source activate" and "source deactivate" by "conda activate" and "conda deactivate" respectively 

### Fixed
- Updated Docker tags in tests

## [0.5.1] - 2022-03-21
### Fixed
- Updated the expected test output files
- Changed the minimum number of samples in the filter_ref_sites_vcf.py test to 1

## [0.5.0] - 2021-01-08
### Changed
- merge_vcf_files.py now sets genotypes with a DHFFC of Inf, Nan, or higher than 4 to missing
- vcf_to_table.py now writes both samples and genotypes into the SAMPLES-VAR column, as opposed to samples only

### Added
- Option to merge_vcf_files.py to merge VCF files containing multiple samples per file, 
- Option to merge_vcf_files.py to specify read depth cutoffs for calling deletions and duplications
- Option to merge_vcf_files.py to generate output that can be processed by [svtyper](https://github.com/hall-lab/svtyper)
- Nextflow scripts to convert output BEDPEs to VCFs and add read depth information to each call using duphold (script heavily inspired by [methylseq](https://github.com/nf-core/methylseq))
- Potentially missing hexdump dependency to Dockerfile

## [0.4.0] - 2020-04-14
### Added
- Option in merge_vcf_files.py to merge VCF files of different samples without genotyping based on read depth computed by duphold
- Option in vcf_to_table.py to calculate number of deletions that are supported by a change in read depth compared to 1000 bp flanking regions
- Option in vcf_to_table.py to generate, for each variant, a list of identifiers of samples that had a non-reference call
- Option in vcf_to_table.py to generate ID field of each variant

### Changed
- Specified conda channels in requirements.yml files to improve reproducibility

### Fixed
- FILTER field is now correctly processed in vcf_to_table.py
- Nextflow binary now has read permissions in the Docker image, so that it can be used by non-root users

## [0.3.0] - 2019-12-16
### Added
- Script to easily filter Hecaton output (nextflow/hecaton_filter_size_cutoff.nf)

### Changed
- Added docker pull option in documentation
- Corrected link to old repository in documentation

### Fixed
- Bug during merging in which homozygous insertions are reported for a sample, even if its VCF gives a reference or no call at that position   
- Typo in hecaton.nf and hecaton_no_align.nf which caused speedseq -m parameter to be equal to the number of threads
- Typos in manual

## [0.2.2] - 2019-09-20
### Added
- Test to check if running Hecaton with pre-generated alignments works correctly

### Changed
- Included additional dependencies in the documentation that are needed when running Hecaton on a fresh Ubuntu install

### Fixed
- Bugs which prevented Hecaton from running properly when using pre-generated alignments

## [0.2.1] - 2019-09-05
### Added
- Quick and dirty local installation script

### Changed
- Made it clearer in the documentation that you need to have installed Nextflow before testing the Hecaton installation
- Updated the version of Hecaton in the Dockerfile, Docker container building script, and Docker installation test script

### Fixed
- Conda environments are now activated properly when running Hecaton with newer versions of Nextflow

## [0.2.0] - 2019-08-17
### Added
- Hecaton can now be run with pre-generated BWA mem alignments
- Certain chromosomal regions can now be excluded from CNV calling
- Script which can convert BEDPE to VCF file (all calls are treated as homozygous)
- Script which merges VCFs of different samples
- Example of nextflow config file for Slurm

### Fixed
- Missing parameter in functional test script

## [0.1.0] - 2019-03-16
### Added
- Initial release of Hecaton

[Unreleased]: https://git.wur.nl/bioinformatics/hecaton/compare/v0.3.0...master
[0.3.0]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.3.0
[0.2.2]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.2.2
[0.2.1]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.2.1
[0.2.0]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.2.0
[0.1.0]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.1.0

