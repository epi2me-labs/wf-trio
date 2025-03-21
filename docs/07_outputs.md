Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Per file read stats | ./fastq_ingress_results/reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Per read stats | ./fastq_ingress_results/reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | ./fastq_ingress_results/reads/fastcat_stats/run_ids | List of run ID's present in reads. | aggregated |
| Meta map json | ./fastq_ingress_results/reads/metamap.json | Meta data used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | ./fastq_ingress_results/reads/{{ alias }}.fastq.gz | Per sample reads concatenated in to one fastq file. | per-sample |
| Trio small variants report | {{ family_id }}.wf-trio-snp-report.html | Analysis and stats relating to the multi-sample small variants VCF. | aggregated |
| Trio SV report | {{ family_id }}.wf-trio-sv-report.html | Analysis and stats relating to the multi-sample SV VCF. | aggregated |
| Individual SV report | {{ sample_id }}.wf-trio-sv-report.html | Analysis and stats relating to the SVs for each individual. | per-sample |
| Individual small variants report | {{ sample_id }}.wf-trio-snp-report.html | Analysis and stats relating to the small variants for each individual. | per-sample |
| Individual QC report | {{ sample_id }}.wf-trio-qc-report.html | Quality control read stats including read length and quality scores. | per-sample |
