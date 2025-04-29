Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Small variant VCF | {{ alias }}.wf_trio_snp.vcf.gz | VCF file with small variant calls (SNPs and indels) for the sample. These will be unphased. | per-sample |
| Small variant GVCF | {{ alias }}.wf_trio_snp.gvcf.gz | GVCF file with the SNPs for the sample. These will be unphased. | per-sample |
| Joint Small variant VCF | {{ family_id }}.wf_trio_snp.vcf.gz | Joint VCF file with Small variant calls (SNPs and indels). | aggregated |
| Structural variant VCF | {{ alias }}.wf_trio_sv.vcf.gz | VCF file with the SVs for the sample. These will be phased if SNP is enabled. | per-sample |
| Joint structural variant VCF | {{ family_id }}.wf_trio_sv.vcf.gz | VCF file with the SVs for all the samples. | aggregated |
| Structural variant SNF | {{ alias }}.wf_trio_sv.snf | SNF file with the SVs for the sample, for onward SV merging. | per-sample |
| Haplotagged alignment file | {{ alias }}.haplotagged.bam | BAM file with the haplotagged reads for the sample. Will be pedigree phased if phased option is used. | per-sample |
| Haplotagged alignment file index | {{ alias }}.haplotagged.bam.bai | The index of the resulting BAM file with the haplotagged reads for the sample. | per-sample |
| RTG Annotation SNP | {{ family_id }}.wf_trio_snp.RTGannot.txt | RTG mendelian inheritance report from SNP VCF. | aggregated |
| RTG Annotation SV | {{ family_id }}.wf_trio_sv.RTGannot.txt | RTG mendelian inheritance report from SV VCF. | aggregated |
| Pedigree file | {{ family_id }}.wf-trio.ped | The pedigree file filtered by family_id. | aggregated |
| Trio small variants report | {{ family_id }}.wf-trio-snp-report.html | Analysis and stats relating to the multi-sample small variants VCF. | aggregated |
| Trio SV report | {{ family_id }}.wf-trio-sv-report.html | Analysis and stats relating to the multi-sample SV VCF. | aggregated |
| Individual SV report | {{ sample_id }}.wf-trio-sv-report.html | Analysis and stats relating to the SVs for each individual. | per-sample |
| Individual small variants report | {{ sample_id }}.wf-trio-snp-report.html | Analysis and stats relating to the small variants for each individual. | per-sample |
| Individual QC report | {{ sample_id }}.wf-trio-qc-report.html | Quality control read stats including read length and quality scores. | per-sample |
