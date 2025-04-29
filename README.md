# Trio Sequencing workflow

Nextflow workflow for trio sequencing analysis.



## Introduction

<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

This repository contains a [nextflow](https://www.nextflow.io/) workflow for analysing variation using trio human genomic data. Specifically the workflow currently performs:

* Small variant calling
* Small variant merging and joint genotyping
* Structural variant calling
* Structural variant merging
* Pedigree phasing



## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 16
+ Memory = 32GB

Approximate run time: Variable, depending on coverage and the individual analyses requested. For 30X coverage, and analysis of phased SNPs and SVs, the workflow will take approximately 12h with recommended resources.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop)
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
	--bed 'wf-trio-demo/Ashkenazi_test.bed' \
	--family_id 'Ashkenazi' \
	--proband_bam 'wf-trio-demo/Ashkenazi_test_data/nf_bam_test_hg002/hg002.bam' \
	--proband_sample_name 'hg002' \
	--pat_bam 'wf-trio-demo/Ashkenazi_test_data/nf_bam_test_hg003/hg003.bam' \
	--pat_sample_name 'hg003' \
	--mat_bam 'wf-trio-demo/Ashkenazi_test_data/nf_bam_test_hg004/hg004.bam' \
	--mat_sample_name 'hg004' \
	--pedigree_file 'wf-trio-demo/ped_file.ped' \
	--ref 'wf-trio-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
	--snp --sv --phased \
	--out_dir demo_results \
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
The `--proband_bam`, `--pat_bam` and `--mat_bam` input parameters for this workflow accept a path to a single BAM file, or a folder containing multiple BAM files for the sample.
Sample names provided to the workflow must match with their corresponding entries in the pedigree file. These names must be specified using the sample_name options: `proband_sample_name`, `pat_sample_name`, `mat_sample_name`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```



## Input parameters

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






## Outputs

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




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
The workflow is composed of 2 subworkflows, each enabled by a command line option:

* [Small variant calling](#3-small-variant-calling-and-joint-genotyping): `--snp`
* [Structural variant calling](#5-structural-variant-calling-and-merging): `--sv`

### 1. Input and data preparation

The workflow requires:
1. A reference genome in [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. Sequencing data for each of the trio samples in the form of a single [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) or a folder of BAM files, either aligned or unaligned. If unaligned input is provided, alignment will be performed using [Minimap2](https://github.com/lh3/minimap2).
3. A [pedigree file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) to describe familial relationships between input samples. Sample names provided to the workflow should match exactly with their entries in the pedigree file. These names must be specified using the sample_name options: `proband_sample_name`, `pat_sample_name`, `mat_sample_name`. Note: the BAM header will be modified by the workflow so the RG ID matches the sample names provided with these options.

The workflow currently only supports the hg38/GRCh38 reference genome. When analysing human data, we recommend using [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz](https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). For more information see [this blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) which outlines potential pitfalls with the various flavours of human references.

The input BAM file can be generated using our [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling/) workflow, which is up to date with the current Dorado releases and models.

### 2. Data QC and pre-processing

The workflow starts by performing multiple checks of the input BAM file, as well as computing:
1. depth of sequencing with [mosdepth](https://github.com/brentp/mosdepth);
2. read alignment statistics with [fastcat](https://github.com/epi2me-labs/fastcat).

Results of these are available in the final report.

### 3. Small variant calling and joint genotyping.

The workflow implements a modified version of [Clair3-Nova](https://github.com/HKU-BAL/Clair3-Nova) to call germline variants in trio.
Clair3-Nova has been modified to improve parallelism in high-performance compute environments.
The workflow will select an appropriate Clair3 and Clair3-Nova model by detecting the basecall model from the input data.
If the input data does not have the required information to determine the basecall model, the workflow will require the basecall model to be provided explicitly with the `--override_basecaller_cfg` option. This will be used by the workflow to automatically determine the Clair3 and Clair3-Nova models to use for small variant calling.

The gVCFs generated by Clair3-Nova are then input in to [GLnexus](https://github.com/dnanexus-rnd/GLnexus) for joint genotyping.

Records in the VCF that are within a homopolymer region greater than 9 base pairs (`HP9`) will be annotated in both the individual and joint small variant VCF using a default BED file internal to the workflow.

### 4. Pedigree phasing 

The workflow can perform pedigree phasing of the variants by using the `--phased` option. The workflow uses [whatshap](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-pedigrees) to perform pedigree phasing of the GLnexus joint variants VCF. Haplotype alleles of a proband are given as paternal|maternal; i.e. the first allele is the one inherited from the father and the second one the allele inherited from the mother. The BAMs will then be haplotagged using the phased VCFs.

### 5. Structural variant calling and merging

The workflow allows for merging of SVs using long-read sequencing data with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).
The workflow will perform SV calling on each individual and the results of these will be used in structural variant merging.
The SV workflow uses an appropriate tandem repeat annotation BED file to improve calling in repetitive regions. You can override the BED file used with the `--tr_bed` parameter. See the [Sniffles](https://github.com/fritzsedlazeck/Sniffles) documentation for more information.

Structural variants can be phased using `--phased`. However, this will cause the workflow to run small variant calling analysis, as SV phasing relies on the haplotagged reads generated in this stage.

Records in the VCF that are within a Tandem repeat (`TR`) region will be annotated in both the individual and joint structural variant VCF using a default BED file internal to the workflow.

### 6. Mendelian inheritance

[RTG Mendelian](https://realtimegenomics.github.io/rtg-tools/rtg_command_reference.html#mendelian) is used to report variants which do not follow Mendelian inheritance. This is run for small variant calling and structural variant multi-sample VCF.




## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-trio/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 

+ *The number of SNVs and indels in the report do not sum up to the number of records, is that normal?* - Yes; this can be due to some multiallelic sites carrying a mixture of SNV and indel alleles.



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




