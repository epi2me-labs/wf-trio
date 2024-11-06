# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Removed
- `--basecaller_cfg` is now redundant as the workflow automatically detects the basecall model from the input data.
- `--clair3_trio_model` is now redundant as the workflow automatically selects the appropriate Clair3 Nova model from the detected basecaller configuration.
### Added
- `--override_basecaller_cfg` parameter allows users to provide a basecall configuration name in cases where automatic basecall model detection fails.
- Error message if more than one basecaller is found across the samples.
- Structural variant calling in trio using `Sniffles2`.
- `--sv` parameter to enable structural variant calling (default: false).
- `--snp` parameter to enable SNP calling with Clair3-Nova (default: false).
### Changed
- Update clair3-nova v0.3.0 to v0.3.1 
- Reconcile template with v5.3.1

## [v0.0.2]
### Fixed
* "failed to chown" errors when pulling Docker images with Singularity.

## [v0.0.1]
### Added
* Workflow provides de novo small variant calling in trio using Clair3-nova.
* Joint variant calling using Glnexus.
* Use RTG tools to checks a multi-sample VCF file for variant calls which do not follow Mendelian inheritance, and compute aggregate sample concordance.
