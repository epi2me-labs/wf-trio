# Trio Sequencing workflow

Nextflow workflow for trio sequencing analysis.



## Introduction

<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

This workflow is not intended to be used by end users.

This workflow can be used for the following:

+ As a template using gitlabs create project from template.
+ For testing of any scripts that are the same across workflows such as scripts in the lib directory.



## Compute requirements

Recommended requirements:

+ CPUs = 2
+ Memory = 2GB

Minimum requirements:

+ CPUs = 2
+ Memory = 2GB

Approximate run time: 5 minutes per sample

ARM processor support: True




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-trio --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-trio
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-trio/wf-trio-demo.tar.gz
tar -xzvf wf-trio-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-trio \
	--proband_bam wf-trio-demo/Ashkenazi_test_TP_data/nf_bam_test_hg002/hg002.bam \
	--pat_bam wf-trio-demo/Ashkenazi_test_TP_data/nf_bam_test_hg003 \
	--mat_bam wf-trio-demo/Ashkenazi_test_TP_data/nf_bam_test_hg004 \
	--ref wf-trio-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| proband_bam | string | Proband aligned BAM to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files for a single sample. |  |
| pat_bam | string | Paternal aligned BAM files to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files for a single sample. |  |
| mat_bam | string | Maternal aligned BAM files to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files for a single sample. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| pedigree_file | string | TSV file containing the columns `family_id`, `individual_id`, `paternal_id`, `maternal_id`, `sex`, `phenotype` | See https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format for further information. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| watch_path | boolean | Enable to continuously watch the input directory for new input files. | This option enables the use of Nextflow’s directory watching feature to constantly monitor input directories for new files. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| proband_sample_name | string | A single sample name for non-multiplexed data. Permissible if passing a single .bam file or directory of .bam files. |  | proband |
| pat_sample_name | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  | paternal |
| mat_sample_name | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  | maternal |
| family_id | string | A single family ID relating to the group of samples. |  | family1 |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| phased | boolean | Perform phasing. | This option enables phasing of SV, SNP and modifications, depending on which sub-workflow has been chosen; see [README](README.md#9-phasing-variants) for more details. | False |
| bed | string | An optional BED file enumerating regions to process for variant calling. |  |  |
| basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | Required for small variant calling. The basecaller configuration is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. | dna_r10.4.1_e8.2_400bps_sup@v4.3.0 |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small and structural variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| GVCF | boolean | Enable to output a gVCF file in addition to the VCF outputs (experimental). | By default the workflow outputs a VCF file containing only records where a variant has been detected. Enabling this option will additionally output a gVCF with records spanning all reference positions regardless of whether a variant was detected in the sample. | True |
| phase_trio | boolean | Output a phased VCF and haplotagged BAM. | By default the workflow will output a phased VCF and haplotagged BAM, set to False to skip. | True |
| glnexus_config | string | glnexus config yaml file, a default one is provided. | See https://github.com/dnanexus-rnd/GLnexus/wiki/Configuration for further details. |  |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-trio-report.html | Report for all samples | aggregated |
| Per file read stats | ./fastq_ingress_results/reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Per read stats | ./fastq_ingress_results/reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | ./fastq_ingress_results/reads/fastcat_stats/run_ids | List of run ID's present in reads. | aggregated |
| Meta map json | ./fastq_ingress_results/reads/metamap.json | Meta data used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | ./fastq_ingress_results/reads/{{ alias }}.fastq.gz | Per sample reads concatenated in to one fastq file. | per-sample |




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.



## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



