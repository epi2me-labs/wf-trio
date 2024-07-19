// Note this is not exactly the same as wf-human-variation.nf (v2.1.0) only required processes are here - all the full variants process is not 
// Throughout as we are using the clair3_nova container we call the clair3 scripts with `python \$CLAIR3_NOVA_PATH/clair3.py
// The 4 other differences are commented with DIFFERENT

import groovy.json.JsonBuilder

def phaser_memory = params.use_longphase ? [8.GB, 32.GB, 56.GB] : [4.GB, 8.GB, 12.GB]
def haptag_memory = [4.GB, 8.GB, 12.GB]

process make_chunks {
    // Do some preliminaries. Ordinarily this would setup a working directory
    // that all other commands would make use off, but all we need here are the
    // list of contigs and chunks.
    label "clair3nova"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta), path(xam), path(xam_idx)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path bed
        path model_path
        val chromosome_codes
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp/CHUNK_LIST", emit: chunks_file
        path "clair_output/tmp/CMD", emit: cmd_file
        path "clair_output/tmp/split_beds", emit: split_beds, optional: true
    script:
        def bedargs = bed.name != 'OPTIONAL_FILE' ? "--bed_fn ${bed}" : ''
        def bedprnt = bed.name != 'OPTIONAL_FILE' ? "--bed_fn=${bed}" : ''
        String ctgs = chromosome_codes.join(',')
        // Define contigs in order to enforce the mitochondrial genome calling, which is otherwise skipped.
        def ctg_name = "--ctg_name ${ctgs}"
        // If a single contig is required, then set it as option
        if (params.ctg_name){
            ctg_name = "--ctg_name ${params.ctg_name}"
        }
        // If all contigs are required, then set the ctg_name to EMPTY
        if (params.include_all_ctgs || bed.name != 'OPTIONAL_FILE'){
            ctg_name = '--ctg_name "EMPTY"'
        }
        // DIFFERENT: removed these as not present in clair3_trio version
        //  --enable_long_indel False \
        //  --cmd_fn clair_output/tmp/CMD
        """
        # CW-2456: save command line to add to VCF file (very long command...)
        mkdir -p clair_output/tmp
        echo "run_clair3.sh --bam_fn=${xam} ${bedprnt} --ref_fn=${ref} --vcf_fn=${params.vcf_fn} --output=clair_output --platform=ont --sample_name=${meta.alias} --model_path=${model_path.simpleName} --ctg_name=${params.ctg_name} ${ctg_name} --include_all_ctgs=${params.include_all_ctgs} --chunk_num=0 --chunk_size=5000000 --qual=${params.min_qual} --var_pct_full=${params.var_pct_full} --ref_pct_full=${params.ref_pct_full} --snp_min_af=${params.snp_min_af} --indel_min_af=${params.indel_min_af} --min_contig_size=${params.min_contig_size}" > clair_output/tmp/CMD
        # CW-2456: prepare other inputs normally
        python \$CLAIR3_NOVA_PATH/clair3.py CheckEnvs \
            --bam_fn ${xam} \
            ${bedargs} \
            --output_fn_prefix clair_output \
            --ref_fn ${ref} \
            --vcf_fn ${params.vcf_fn} \
            ${ctg_name} \
            --chunk_num 0 \
            --chunk_size 5000000 \
            --include_all_ctgs ${params.include_all_ctgs} \
            --threads 1  \
            --qual ${params.min_qual} \
            --sampleName ${meta.alias} \
            --var_pct_full ${params.var_pct_full} \
            --ref_pct_full ${params.ref_pct_full} \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} 
        """
}


process pileup_variants {
    // Calls variants per region ("chunk") using pileup network.
    label "clair3nova"
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy = {task.exitStatus in [134,137,140] ? 'retry' : 'finish'}
    maxRetries 1
    input:
        each region
        tuple val(meta), path(xam), path(xam_idx)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path model
        path bed
        path command
        path "split_bed"
    output:
        // TODO: make this explicit, why is pileup VCF optional?
        path "pileup_*.vcf", optional: true, emit: pileup_vcf_chunks
        path "gvcf_tmp_path/*", optional: true, emit: pileup_gvcf_chunks
    script:
        // note: the VCF output here is required to use the contig
        //       name since that's parsed in the SortVcf step
        // note: snp_min_af and indel_min_af have an impact on performance
        def bedargs = bed.name != 'OPTIONAL_FILE' ? "--bed_fn ${bed} --extend_bed split_bed/${region.contig}" : ''
        //  DIFFERENT: removed --cmd_fn ${command} \
        // Also we use callvarbam instead of callvariantsfromcffi
        // if no new line at end of vcf added one as sometimes it was truncated and caused mangled vcf on concat
        """
        python \$CLAIR3_NOVA_PATH/clair3.py CallVarBam \
            --chkpnt_fn ${model}/pileup \
            --bam_fn ${xam} \
            --call_fn pileup_${region.contig}_${region.chunk_id}.vcf \
            --ref_fn ${ref} \
            --ctgName ${region.contig} \
            --chunk_id ${region.chunk_id} \
            --chunk_num ${region.total_chunks} \
            --platform ont \
            --sampleName ${meta.alias} \
            --fast_mode False \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --call_snp_only False \
            --gvcf True \
            --enable_long_indel False \
            --temp_file_dir gvcf_tmp_path \
            --pileup \
            ${bedargs} \
            --base_err ${params.base_err} \
            --gq_bin_size ${params.gq_bin_size} \
            --keep_iupac_bases False
        if [ -f "pileup_${region.contig}_${region.chunk_id}.vcf" ]; then
            sed -i -e '\$a\\' "pileup_${region.contig}_${region.chunk_id}.vcf"
        fi
        
        """
}


process aggregate_pileup_variants {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network. Determines quality filter for selecting variants
    // to use for phasing.
    // DIFFERENT: added --var_pct_phasing to select qual
    label "clair3nova"
    cpus 2
    memory 4.GB
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        // these need to be named as original, as program uses info from
        // contigs file to filter
        path "input_vcfs/*"
        path contigs
        path command
        val meta
    output:
        tuple val(meta), path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi"), emit: pileup_vcf
        path "phase_qual", emit: phase_qual
    shell:
        '''
        pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf \
            --input_dir input_vcfs/ \
            --vcf_fn_prefix pileup \
            --output_fn pileup.vcf \
            --sampleName !{meta.alias} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        # Replaced bgzip with the faster bcftools index -n
        if [ "$( bcftools index -n pileup.vcf.gz )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi

        bgzip -@ !{task.cpus} -fdc pileup.vcf.gz | \
            pypy \$CLAIR3_NOVA_PATH/clair3.py SelectQual --phase --var_pct_phasing !{params.var_pct_phasing} --output_fn .
        '''
}


process select_het_snps {
    // Filters a VCF by contig, selecting only het SNPs.
    label "clair3nova"
    cpus 2
    memory 4.GB
    input:
        each contig
        tuple val(meta), path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectHetSnp.py#L29
        path "split_folder/phase_qual"
    output:
        tuple val(contig), path("split_folder/${contig}.vcf.gz"), path("split_folder/${contig}.vcf.gz.tbi"), emit: het_snps_vcf
    shell:
        '''
        pypy \$CLAIR3_NOVA_PATH/clair3.py SelectHetSnp \
            --vcf_fn pileup.vcf.gz \
            --split_folder split_folder \
            --ctgName !{contig}

        bgzip -c split_folder/!{contig}.vcf > split_folder/!{contig}.vcf.gz
        tabix split_folder/!{contig}.vcf.gz
        '''
}

process phase_contig {
    // Tags reads in an input BAM from heterozygous SNPs
    // The haplotag step was removed in clair-v0.1.11 so this step re-emits
    //   the original BAM and BAI as phased_bam for compatability,
    //   but adds the VCF as it is now tagged with phasing information
    //   used later in the full-alignment model
    cpus 4
    memory { phaser_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    label "clair3nova"
    input:
        tuple val(contig), path("snps.gz"), path("snps.gz.tbi"), val(meta), path(xam), path(xam_idx), path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple val(contig), path(xam), path(xam_idx), path("phased_${contig}.vcf.gz"), path("phased_${contig}.vcf.gz.tbi"), emit: phased_bam_and_vcf
    script:
        def output_threads = Math.max(task.cpus - 1, 1)
        if (params.use_longphase)
        """
        echo "Using longphase for phasing"
        longphase phase --ont -o phased_${contig} \
            -s ${het_snps} -b ${xam} -r ${ref} -t ${task.cpus}
        bgzip phased_${contig}.vcf
        tabix -f -p vcf phased_${contig}.vcf.gz
        """
        else
        """
        echo "Using whatshap for phasing"
       
        whatshap phase \
            --output phased_${contig}.vcf.gz \
            --reference ${ref} \
            --chromosome ${contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            "snps.gz" \
            ${xam}
        tabix -f -p vcf phased_${contig}.vcf.gz
     
        """
}


process post_clair_contig_haplotag {
    // Tags reads in an input BAM from heterozygous SNPs
    // Also haplotag for those modes that need it
    // We emit BAM as the STR workflow does not fully support CRAM, and so the
    // STR workflow can start while the haplotagged XAM is being catted and
    // written for the final output
    label "clair3nova"
    cpus 4
    // Define memory from phasing tool and number of attempt
    memory { haptag_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple val(contig),
            path(xam), path(xam_idx),
            path(vcf), path(tbi),
            path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
            val(meta)
    output:
        tuple val(contig), val(meta), path("${contig}_hp.bam"), path("${contig}_hp.bam.bai"), emit: phased_bam
    script:
    """
    whatshap haplotag \
        --reference ${ref} \
        --ignore-read-groups \
        --regions ${contig} \
        phased_${contig}.vcf.gz \
        ${xam} \
    | samtools view -O bam --reference $ref -@3 -o ${contig}_hp.bam##idx##${contig}_hp.bam.bai --write-index
    """
}

// DIFFERENT: Additional process to aggregate the pileup gvcf required by callvarbam_trio
process aggregate_pileup_gvcf {
    label "clair3nova"
    maxRetries 3
    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path "merge_outputs_gvcf/*"
        path contigs
        val(meta)
    output:
        tuple val(meta), path("non_var.gvcf"), emit: agg_gvcf
    script:
        """
        pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf \
        --input_dir merge_outputs_gvcf \
        --vcf_fn_suffix ".tmp.gvcf" \
        --sampleName ${meta.alias} \
        --output_fn non_var.gvcf \
        --ref_fn ${ref} \
        --contigs_fn ${contigs}
        """
}


// This is a hilarious trick from CW to present models inside the container as
// outside the container, by exporting them out of the container back to workdir.
// This saves us passing around tuples of val(inside) and path(outside).
process lookup_clair3_model {
    label "clair3nova"
    cpus 1
    memory 2.GB
    input:
        path("lookup_table")
        val basecall_model
    output:
        tuple env(clair3_model), path("model/")
    shell:
    '''
    clair3_model=$(resolve_clair3_model.py lookup_table '!{basecall_model}' clair3)
    cp -r ${CLAIR_MODELS_PATH}/${clair3_model} model
    echo "Basecall model: !{basecall_model}"
    echo "Clair3 model  : ${clair3_model}"
    '''
}


process lookup_clair3_nova_model {
    label "clair3nova"
    cpus 1
    memory 2.GB
    input:
        path("lookup_table")
        val basecall_model
    output:
        tuple env(clair3_nova_model), path("nova_model/")
    shell:
    '''
    clair3_nova_model=$(resolve_clair3_model.py lookup_table '!{basecall_model}' clair3_nova)
    cp -r ${CLAIR_NOVA_MODELS_PATH}/${clair3_nova_model} nova_model
    echo "Basecall model: !{basecall_model}"
    echo "Clair3 nova model  : ${clair3_nova_model}"
    '''
}