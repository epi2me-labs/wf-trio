# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.9.1]
This patch release of wf-trio addresses a problem identified in Sniffles 2.6.2, which was introduced by our wf-trio 0.9.0 release. All users of wf-trio 0.9.0 should immediately adopt this release.

We have identified that in some cases Sniffles 2.6.2 consumes much more RAM than it has been allocated, and is killed by the system without raising an error, resulting in an empty or possibly truncated structural variant VCF. We have worked on [a patch with the Sniffles developers](https://github.com/fritzsedlazeck/Sniffles/pull/565) to ensure that the Sniffles program will exit with a fatal error in this scenario; allowing the workflow to safely terminate without providing an empty or truncated VCF.

This patch is provided to ensure that in this case, users do not unknowingly analyse empty or truncated VCFs, but users should expect a further update to wf-trio to address the underlying Sniffles memory problem.

### Changed
- Sniffles 2.6.2 patched to incorporate [Sniffles/#565](https://github.com/fritzsedlazeck/Sniffles/pull/565) to avoid empty SV VCF if all analysis workers are killed by the operating system.
### Removed
- wf-trio 0.9.0 has been untagged to prevent pipelines pinning this version.

## [v0.9.0]
wf-trio is a new workflow for the analysis of combined whole genome sequencing data from proband and both parents.
This open release supersedes our closed beta and is focused on the calling of small variants, structural variants, and pedigree phasing of variants for the discovery and characterisation of proband variants.
We look forward to seeing how our workflow supports your research.