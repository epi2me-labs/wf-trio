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
- Handling for v4.3.0 HAC and v5.0.0 HAC and SUP basecaller configurations.
- Include the pedigree file in the outputs.
- Alignment resource parameters `ubam_map_threads`, `ubam_sort_threads`, `ubam_bam2fq_threads`.
### Changed
- Phasing performed is now pedigree phasing.
- Update whatshap v1.7.0 to v2.3.0. 
- Update clair3-nova v0.3.0 to v0.3.2 which reduces incomplete variant calls and adds the reference calls to GVCF outputs.
- Transfer De novo probabilities from INFO to FORMAT field.
- Modified the phased Joint SNP VCF header to reflect replacing `Number=.` with `Number=R` in the AD field.
- `{alias}.wf_trio_snp.vcf.gz` output files are now the VCF's output by clair3-nova.
- `{family_id}.wf_trio_snp.vcf.gz` output file contains the joint pedigree phased variant calling.
- A default tandem repeat `tr_bed` is provided to sniffles2. Can be overwritten with the `tr_bed` parameter.
- Annotation of tandem repeat regions and long homopolymers with `bcftools annotate` now uses ID, REF and ALT column to avoid records with partial overlaps to be annotated.
- Output per sample files in per sample folders.
- Updated to wf-template v5.6.1, changing:
    - Reduce verbosity of debug logging from fastcat which can occasionally occlude errors found in FASTQ files during ingress.
    - Log banner art to say "EPI2ME" instead of "EPI2ME Labs" to match current branding. This has no effect on the workflow outputs.
    - pre-commit configuration to resolve an internal dependency problem with flake8. This has no effect on the workflow.
### Fixed
- Updated to wf-template v5.6.1, fixing:
    - dacite.exceptions.WrongTypeError during report generation when barcode is null.
    - Sequence summary read length N50 incorrectly displayed minimum read length, it now correctly shows the N50.
    - Sequence summary component alignment and coverage plots failed to plot under some conditions.


## [v0.0.2]
### Fixed
* "failed to chown" errors when pulling Docker images with Singularity.

## [v0.0.1]
### Added
* Workflow provides de novo small variant calling in trio using Clair3-nova.
* Joint variant calling using Glnexus.
* Use RTG tools to checks a multi-sample VCF file for variant calls which do not follow Mendelian inheritance, and compute aggregate sample concordance.
