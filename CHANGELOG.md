# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed
- Added docker pull option in documentation

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

[Unreleased]: https://git.wur.nl/bioinformatics/hecaton/compare/v0.2.1...master
[0.2.1]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.2.1
[0.2.0]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.2.0
[0.1.0]: https://git.wur.nl/bioinformatics/hecaton/tags/v0.1.0

