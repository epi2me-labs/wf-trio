// Trio

def phaser_memory = params.use_longphase ? [8.GB, 32.GB, 56.GB] : [4.GB, 8.GB, 12.GB]

process checkEnv_trio {
    label "clair3nova"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta_proband), path("proband.bam"), path("proband.bam.bai")
        tuple val(meta_pat), path("pat.bam"), path("pat.bam.bai")
        tuple val(meta_mat), path("mat.bam"), path("mat.bam.bai")
        path bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val chromosome_codes
        path model_path
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp", emit: tmp
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
    """
    mkdir -p clair_output/tmp
    echo "run_clair3.sh \
        --bam_fn="proband.bam" \
        ${bedprnt} \
        --ref_fn=${ref} \
        --vcf_fn=${params.vcf_fn} \
        --output=clair_output \
        --platform=ont \
        --sample_name=${meta_proband.alias} \
        --model_path=${model_path.simpleName} \
        ${ctg_name} \
        --chunk_num=0 \
        --chunk_size=5000000 \
        --qual=${params.min_qual} \
        --var_pct_full=${params.var_pct_full} \
        --ref_pct_full=${params.ref_pct_full} \
        --snp_min_af=${params.snp_min_af} \
        --indel_min_af=${params.indel_min_af} \
        --min_contig_size=${params.min_contig_size} \
        --is_denovo False" > clair_output/tmp/CMD

    python \$CLAIR3_NOVA_PATH/clair3.py CheckEnvs_Trio \
        --bam_fn_c proband.bam \
        --bam_fn_p1 mat.bam \
        --bam_fn_p2 pat.bam \
        --output_fn_prefix clair_output \
        ${bedargs} \
        --ref_fn ${ref} \
        --vcf_fn ${params.vcf_fn} \
        ${ctg_name} \
        --chunk_num 0 \
        --chunk_size 5000000 \
        --include_all_ctgs  ${params.include_all_ctgs} \
        --threads 1 \
        --qual ${params.min_qual} \
        --sampleName_c  ${meta_proband.alias} \
        --sampleName_p1 ${meta_pat.alias}  \
        --sampleName_p2 ${meta_mat.alias}  \
        --var_pct_full ${params.var_pct_full} \
        --ref_pct_full ${params.ref_pct_full} \
        --snp_min_af ${params.snp_min_af} \
        --indel_min_af ${params.indel_min_af}
    """
}


process selectCandidates_trio {
    label "clair3nova"
    cpus 2
    memory 4.GB
    input:
        each contig
        tuple val(meta_proband), path("proband_sample.vcf.gz"), path("proband_sample.vcf.gz.tbi")
        tuple val(meta_pat), path("pat_sample.vcf.gz"), path("pat_sample.vcf.gz.tbi")
        tuple val(meta_mat), path("mat_sample.vcf.gz"), path("mat_sample.vcf.gz.tbi")
    output:
       tuple val(contig), path("candidate_bed/${contig}.*"), emit: candidate_bed
       tuple val(contig), path("candidate_bed"), emit: candidate_bed_folder
    script:
    // Question for Clair3 devs: Why are pct full and var pctfull hard coded in the trio script? Is this a mistake?
    """
    mkdir candidate_bed
    pypy \$CLAIR3_NOVA_PATH/clair3.py SelectCandidates_Trio \
    --alt_fn_c proband_sample.vcf.gz \
    --alt_fn_p1 pat_sample.vcf.gz \
    --alt_fn_p2 mat_sample.vcf.gz \
    --candidate_bed candidate_bed \
    --sampleName trio \
    --ref_pct_full 0.03 \
    --var_pct_full 1.0 \
    --ref_var_max_ratio 1.2 \
    --ctgName ${contig}
    """
}

// Spuriously  pileups can take more than 4.GB of memory and up to 8.GB
// so retry with more memory if it fails
process callVarBam_trio {
    label "clair3nova"
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy = {task.exitStatus in [134,137,140,142] ? 'retry' : 'finish'}
    maxRetries 1
    input:
       tuple val(contig),
        val(meta_proband), path("proband.bam"), path("proband.bam.bai"),
        val(meta_pat), path("pat.bam"), path("pat.bam.bai"),
        val(meta_mat), path("mat.bam"), path("mat.bam.bai"),
        path(candidate_bed),
        path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
        path("tmp"), path(model_path)
    output:
       tuple val(contig), val(meta_proband), path("*${meta_proband.alias}_${contig}*.vcf"), emit: proband_calls
       tuple val(contig), val(meta_pat), path("*${meta_pat.alias}_${contig}*.vcf"), emit: pat_calls
       tuple val(contig), val(meta_mat), path("*${meta_mat.alias}_${contig}*.vcf"), emit: mat_calls
    script:
       def region = candidate_bed.name
    // TODO: Add error if no variants output at this point as currently fails silently
    """
    mkdir -p calls
    echo $contig
    python \$CLAIR3_NOVA_PATH/clair3.py CallVarBam_Denovo \
    --chkpnt_fn "${model_path}/nova" \
    --bam_fn_c proband.bam \
    --bam_fn_p1 pat.bam \
    --bam_fn_p2 mat.bam \
    --sampleName_c  ${meta_proband.alias} \
    --sampleName_p1 ${meta_pat.alias}  \
    --sampleName_p2 ${meta_mat.alias}  \
    --call_fn_c trio_${meta_proband.alias}_${contig}_${region}.vcf \
    --call_fn_p1 trio_${meta_pat.alias}_${contig}_${region}.vcf \
    --call_fn_p2 trio_${meta_mat.alias}_${contig}_${region}.vcf \
    --use_gpu 0 \
    --ref_fn ${ref} \
    --ctgName ${contig} \
    --platform ont \
    --full_aln_regions "${candidate_bed}" \
    --gvcf True \
    --showRef True \
    --tmp_path tmp \
    --phasing_info_in_bam \
    --keep_iupac_bases False
    """
}

// All vcf contigs per sample are merged and sorted
process sortVcf_Trio {
    label "clair3nova"
    cpus 2
    maxRetries 3
    memory { 8.GB * task.attempt }
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
       tuple val(meta), path("calls/*"),
        path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
        path("tmp")
    output:
       tuple val(meta), path("${meta.alias}_c3t.vcf.gz"), path("${meta.alias}_c3t.vcf.gz.tbi"), emit: trio_pileup
    script:
    """
    pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf_Trio \
        --input_dir calls \
        --vcf_fn_prefix  "trio" \
        --output_fn "${meta.alias}_c3t.vcf" \
        --sampleName "${meta.alias}" \
        --ref_fn $ref \
        --tmp_path tmp \
        --contigs_fn tmp/CONTIGS
    """
}


process mergeVcf_Trio {
    label "clair3nova"
    cpus 1
    input:
        tuple val(contig), val(meta),
            path("${meta.alias}_c3t.vcf.gz"),
            path("${meta.alias}_c3t.vcf.gz.tbi"),
            path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi"),
            path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
            path(aggregated_gvcf), path("candidate_beds")
    output:
       tuple val(meta), val(contig), path("merge_${contig}.vcf"), emit: merged_vcf
       tuple val(meta), val(contig), path("merge_${contig}.gvcf"), emit: merged_gvcf
    script:
        """
        pypy \$CLAIR3_NOVA_PATH/clair3.py MergeVcf_Trio \
        --pileup_vcf_fn pileup.vcf.gz \
        --trio_vcf_fn ${meta.alias}_c3t.vcf.gz \
        --output_fn merge_${contig}.vcf \
        --gvcf_fn  "merge_${contig}.gvcf.lz4" \
        --platform ont \
        --print_ref_calls True \
        --gvcf True \
        --bed_fn_prefix candidate_beds \
        --haploid_precise False \
        --haploid_sensitive False \
        --non_var_gvcf_fn ${aggregated_gvcf} \
        --ref_fn ${ref} \
        --ctgName ${contig}

        # mergevcf trio outputs gvcf as lz4 as is an intermediate file within clair3-nova
        # to make it compatible with downstream tools in this workflow mainly glnexus uncompress
        lz4 "merge_${contig}.gvcf.lz4" "merge_${contig}.gvcf"
       
        """
}


process sortmergedVcf_Trio {
    label "clair3nova"
    cpus 2
    maxRetries 3
    memory { 8.GB * task.attempt }
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
       tuple val(meta), path("merged_calls/*"),
        path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
        path("tmp")
    output:
       tuple val(meta), path("${meta.alias}.wf_trio_snp.vcf.gz"), path("${meta.alias}.wf_trio_snp.vcf.gz.tbi"), emit: sorted_trio
    script:
    """
    pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf_Trio \
    --input_dir merged_calls \
    --vcf_fn_prefix  "merge" \
    --output_fn "${meta.alias}.wf_trio_snp.vcf" \
    --sampleName "${meta.alias}" \
    --ref_fn "$ref" \
    --tmp_path tmp \
    --contigs_fn tmp/CONTIGS
    """
}


process sortmergedGVcf_Trio {
    // all gvcf contigs per sample are merged and sorted
    // merge gvcf calls include gvcfs from clair3 pile up
    label "clair3nova"
    cpus 2
    maxRetries 3
    memory { 8.GB * task.attempt }
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
       tuple val(meta), path("merged_gvcf_calls/*"), path(ref), path(ref_idx), path(ref_cache), env(REF_PATH), path("tmp")
    output:
       tuple val(meta), path("${meta.alias}.wf_trio_snp.gvcf.gz"), path("${meta.alias}.wf_trio_snp.gvcf.gz.tbi"), emit: sorted_trio
    script:
    """
    pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf_Trio \
    --input_dir merged_gvcf_calls \
    --vcf_fn_suffix ".gvcf" \
    --output_fn "${meta.alias}.wf_trio_snp.gvcf" \
    --sampleName $meta.alias \
    --ref_fn $ref \
    --tmp_path tmp \
    --contigs_fn tmp/CONTIGS
    """
}


process phase_trio {
    tag "$contig:${meta.alias}"
    label "clair3nova"
    cpus 4
    memory { phaser_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
       tuple val(meta),val(contig),
        path(het_snps),
        path(xam), path(xam_idx),
        path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
       tuple val(meta),val(contig), path("phased_${contig}.vcf"), emit: phased
    script:
    """
    whatshap phase \
        --output phased_${contig}.vcf \
        --reference ${ref} \
        --ignore-read-groups \
            ${het_snps} \
            ${xam}
    """

}


process sortphased_Trio {
    label "clair3nova"
    cpus 2
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    memory { 8.GB * task.attempt }
    input:
       tuple val(meta), path("phased_calls/*"),
        path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
        path("tmp")
    output:
       tuple val(meta), path("${meta.alias}.phased.vcf.gz"), path("${meta.alias}.phased.vcf.gz.tbi"), emit: final_vcf
    script:
    """
    pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf_Trio \
    --input_dir phased_calls \
    --vcf_fn_prefix  "phased" \
    --output_fn "${meta.alias}".phased.vcf \
    --sampleName $meta.alias \
    --ref_fn $ref \
    --tmp_path tmp \
    --contigs_fn tmp/CONTIGS
    """
}


process haplotag_trio {
    cpus 4
    memory 4.GB
    label "clair3nova"
    input:
       tuple val(meta), val(contig), path(het_snps),
        path(xam), path(xam_idx),
        path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
       tuple val(meta), path("${contig}_hp.bam"), emit: haplotag_bam
    script:
    """
    bgzip phased_${contig}.vcf 
    tabix -f -p vcf phased_${contig}.vcf.gz
    whatshap haplotag \
        --reference ${ref} \
        --regions ${contig} \
        --ignore-read-groups \
        phased_${contig}.vcf.gz  \
        ${xam} \
    | samtools view -O bam --reference ${ref} -@3 -o ${contig}_hp.bam##idx##${contig}_hp.bam.bai --write-index
    """
}


process cat_haplotagged_contigs {
    label "wf_human_snp"
    cpus 4
    memory 15.GB // cat should not need this, but weirdness occasionally strikes
    input:
        tuple val(meta), path("contig_bams/*"),
            path(ref), path(ref_idx), path(ref_cache), env(REF_PATH) // intermediate input always BAM here
        tuple val(xam_fmt), val(xai_fmt) // currently always .bam, .bai. But will add .cram, crai in a future release
    output:
        tuple path("${meta.alias}.haplotagged.${xam_fmt}"), path("${meta.alias}.haplotagged.${xam_fmt}.${xai_fmt}"), emit: merged_xam
    script:
    def threads = task.cpus - 1
    """
    # ensure this bit is idempotent as it will inevitably not be called so
    if [ -f seq_list.txt ]; then
        rm seq_list.txt
    fi
    # pick the "first" bam and read its SQ list to determine sort order
    samtools view -H --no-PG `ls contig_bams/*_hp.bam | head -n1` | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*,\\1,p' > seq_list.txt
    # append present contigs to a file of file names, to cat according to SQ order
    while read sq; do
        if [ -f "contig_bams/\${sq}_hp.bam" ]; then
            echo "contig_bams/\${sq}_hp.bam" >> cat.fofn
        fi
    done < seq_list.txt
    if [ ! -s cat.fofn ]; then
        echo "No haplotagged inputs to cat? Are the input file names staged correctly?"
        exit 70 # EX_SOFTWARE
    fi

    # cat just cats, if we want bam, we'll have to deal with that ourselves
    if [ "${xam_fmt}" = "cram" ]; then
        samtools cat -b cat.fofn --no-PG -o - | samtools view --no-PG -@ ${threads} --reference ${ref} -O CRAM --write-index -o "${meta.alias}.haplotagged.cram##idx##${meta.alias}.haplotagged.cram.crai"
    else
        samtools cat -b cat.fofn --no-PG -@ ${threads} -o "${meta.alias}.haplotagged.bam"
        samtools index -@ ${threads} -b "${meta.alias}.haplotagged.bam"
    fi
    """
}


process sortgvcf_Trio_contig {
    label "clair3nova"
    cpus 2
    maxRetries 3
    memory { 8.GB * task.attempt }
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
       tuple val(meta), val(contig), path("calls/*"), path(ref), path(ref_idx), path(ref_cache), env(REF_PATH), path("tmp")
    output:
       tuple val(meta), val(contig), path("${meta.alias}.${contig}.gvcf.gz"), path("${meta.alias}.${contig}.gvcf.gz.tbi")
    script:
    """
    echo -e ${contig} > CONTIGS
    pypy \$CLAIR3_NOVA_PATH/clair3.py SortVcf_Trio \
    --input_dir calls \
    --vcf_fn_suffix ".gvcf" \
    --output_fn "${meta.alias}.${contig}.gvcf" \
    --sampleName "${meta.alias}" \
    --ref_fn $ref \
    --tmp_path tmp \
    --contigs_fn CONTIGS
    """
}


process glnexus {
    cpus 6
    memory { 8.GB * task.attempt }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    label "wftrio"
    input:
       tuple val(contig),
       path("proband.gvcf.gz"), path("paternal.gvcf.gz"), path("maternal.gvcf.gz"),
       path("proband.vcf"), val(family_id), path("gl_config.yml")
    output:
       tuple val(family_id), path("${family_id}.${contig}.vcf.gz"), emit: gl_vcf
    script:
    def threads = task.cpus - 2
    """
    glnexus_cli --threads ${threads} --dir glnexus_DB --config "gl_config.yml" \
    "proband.gvcf.gz" "paternal.gvcf.gz"  "maternal.gvcf.gz"  \
    | bcftools view --threads 1 -o merged.vcf.gz -O z - && bcftools index --threads 1 -t merged.vcf.gz

    echo -e ${params.proband_sample_name}"\\n"${params.pat_sample_name}"\\n"${params.mat_sample_name} > samples.lst
    bgzip "proband.vcf"
    tabix -p vcf "proband.vcf.gz"

    bcftools annotate --collapse none -c "INFO/DNP:=DNP" -a "proband.vcf.gz" merged.vcf.gz \
    | bcftools view -S samples.lst -o "${family_id}.${contig}.vcf.gz" && bcftools index --threads 1 -t "${family_id}.${contig}.vcf.gz"

    """
}


process rtgTools {
    cpus 4
    memory 16.GB
    label "wftrio"
    input:
       tuple val(family_id), path("joint_vcf.vcf.gz"), path("joint_vcf.vcf.gz.tbi")
       path(ref_file)
       path(ped_file)
    output:
       tuple val(family_id), path("*_RTGannot.txt"), emit: summary
    script:
    """
    rtg RTG_MEM=16G format -o reference.sdf $ref_file

    rtg RTG_MEM=15G mendelian -i "joint_vcf.vcf.gz" -o RTGannot.vcf.gz \
     --pedigree=${ped_file} \
     -t reference.sdf | tee ${family_id}_RTGannot.txt

    """
}
