### Workflow Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| snp | boolean | Call small variants | If this option is selected, small variant calling will be carried out using Clair3-Nova, and merging and joint genotyping of variants will be carried out using Glnexus. | True |
| sv | boolean | Call structural variants. | If this option is selected, structural variant calling and merging will be carried out using Sniffles2. | False |
| phased | boolean | This option enables pedigree phasing. | If true the workflow will output the joint phased VCF and haplotagged BAMs. Set to false to skip. Note: Individual SNP VCFs will be unphased. | True |


### Main Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| proband_bam | string | Proband aligned or unaligned BAM to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files for a single sample. |  |
| proband_sample_name | string | A sample name for the proband sample. | Must match the sample name provided in the pedigree input file. |  |
| pat_bam | string | Paternal aligned or unaligned BAM to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files for a single sample. |  |
| pat_sample_name | string | A sample name for the paternal sample. | Must match the sample name provided in the pedigree input file. |  |
| mat_bam | string | Maternal aligned or unaligned BAM to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files for a single sample. |  |
| mat_sample_name | string | A sample name for the maternal sample. | Must match the sample name provided in the pedigree input file. |  |
| family_id | string | A single family ID relating to the group of samples. | Should match the family name provided in the pedigree input file. | family1 |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| pedigree_file | string | A pedigree file describing the familial relationships between input samples in the [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). | Provide a TSV file containing the columns `family_id`, `individual_id`, `paternal_id`, `maternal_id`, `sex`, `phenotype` with relevant values in the rows. Note: If additional records are present in the file only the ones specified by the sample name parameters will be used. |  |
| bed | string | An optional BED file enumerating regions to process for variant calling. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Small variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| include_all_ctgs | boolean | Call variants on all sequences in the reference, otherwise small and structural variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| override_basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data. Providing a value for this option will override the model found in the data. If the model cannot be found in the header, it must be provided with this option. The basecaller model is used to automatically select the appropriate small variant calling model for Clair3 and Clair3-Nova. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. Note: Clair3 v4.3.0 models will be used for v5.0.0 basecalled data.  |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| client_fields | string | A JSON file of key value pairs to display on the report. | This is used to populate a table of additional information (for example, about upstream components of an assay) to the workflow report. |  |


