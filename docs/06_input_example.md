<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
The `--proband_bam`, `--pat_bam` and `--mat_bam` input parameters for this workflow accept a path to a single BAM file, or a folder containing multiple BAM files for the sample.
Sample names provided to the workflow must match with their corresponding entries in the pedigree file. These names must be specified using the sample_name options: `proband_sample_name`, `pat_sample_name`, `mat_sample_name`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```