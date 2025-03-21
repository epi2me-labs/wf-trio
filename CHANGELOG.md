# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Removed
- `--basecaller_cfg` is now redundant as the workflow automatically detects the basecall model from the input data.
- `--clair3_trio_model` is now redundant as the workflow automatically selects the appropriate Clair3 Nova model from the detected basecaller configuration.
- `--phase_trio` has been removed as pedigree phasing is automatically performed when using `--phased`.
### Added
- `--override_basecaller_cfg` parameter allows users to provide a basecall configuration name in cases where automatic basecall model detection fails.
- Error message if more than one basecaller is found across the samples.
- Structural variant calling in trio using `Sniffles2`.
- `--sv` parameter to enable structural variant calling (default: false).
- `--snp` parameter to enable SNP calling with Clair3-Nova (default: false).
- RefCall SNPs are now filtered out of the final SNP VCFs.
- `--phased` parameter to enable pedigree phasing of the SNP VCF's. (default: true) 
- Output additional reports including Joint SNP report, Joint SV report and individual SNP and SV reports.
- Annotation of tandem repeat regions in individual and joint SV VCF's, will use default BED file internal to the workflow.
- Annotation of long homopolymers in individual VCF's/GVCF's and joint SNP VCF's, will use default BED file internal to the workflow.
### Changed
- Phasing performed is now pedigree phasing.
- Update whatshap v1.7.0 to v2.3.0. 
- Update clair3-nova v0.3.0 to v0.3.2 which reduces incomplete variant calls and adds the reference calls to GVCF outputs.
- Transfer De novo probabilities from INFO to FORMAT field.
- Modified the phased Joint SNP VCF header to reflect replacing `Number=.` with `Number=R` in the AD field.
- `{alias}.wf_trio_snp.vcf.gz` output files are now the VCF's output by clair3-nova.
- `{family_id}.wf_trio_snp.vcf.gz` output file contains the joint pedigree phased variant calling.
- A default tandem repeat `tr_bed` is provided to sniffles2. Can be overwritten with the `tr_bed` parameter.
- Reconciled workflow with wf-template v5.5.0.

## [v0.0.2]
### Fixed
* "failed to chown" errors when pulling Docker images with Singularity.

## [v0.0.1]
### Added
* Workflow provides de novo small variant calling in trio using Clair3-nova.
* Joint variant calling using Glnexus.
* Use RTG tools to checks a multi-sample VCF file for variant calls which do not follow Mendelian inheritance, and compute aggregate sample concordance.
