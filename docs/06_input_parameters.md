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
| override_basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data. Providing a value for this option will override the model found in the data. If the model cannot be found in the header, it must be provided with this option. The basecaller model is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. |  |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small and structural variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| GVCF | boolean | Enable to output a gVCF file in addition to the VCF outputs (experimental). | By default the workflow outputs a VCF file containing only records where a variant has been detected. Enabling this option will additionally output a gVCF with records spanning all reference positions regardless of whether a variant was detected in the sample. | True |
| phase_trio | boolean | Output a phased VCF and haplotagged BAM. | By default the workflow will output a phased VCF and haplotagged BAM, set to False to skip. | True |
| glnexus_config | string | glnexus config yaml file, a default one is provided. | See https://github.com/dnanexus-rnd/GLnexus/wiki/Configuration for further details. |  |


