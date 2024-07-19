include { 
    selectCandidates_trio;
    callVarBam_trio;
    sortVcf_Trio;
    mergeVcf_Trio;
    sortmergedVcf_Trio;
    sortmergedGVcf_Trio;
    phase_trio;
    sortphased_Trio;
    haplotag_trio;
    cat_haplotagged_contigs;
    checkEnv_trio;
    glnexus;
    rtgTools;
    sortgvcf_Trio_contig;
} from '../modules/local/wf-trio'

include {
    index_ref_fai;
    cram_cache;
    concat_vcfs;
} from '../modules/local/common'

include {
    snp as snp_proband; snp as snp_mat; snp as snp_pat;
} from "../subworkflows/wf-human-snp"

OPTIONAL_FILE = "$projectDir/data/OPTIONAL_FILE"
// workflow module
workflow trio {
    take:
        proband_bam
        pat_bam
        mat_bam
        snp_bed
        clair3_model
        clair3_nova_model
        fam_id
        gl_conf
        ped_file

    main:

        // To do: Add cram handling and CI test
        extensions = ['bam', 'bai']
        
        // Programmatically define chromosome codes.
        ArrayList chromosome_codes = []
        ArrayList chromosomes = [1..22] + ["X", "Y", "M", "MT"]
        for (N in chromosomes.flatten()){
            chromosome_codes += ["chr${N}", "${N}"]
        }

        // Prepare bams required for snp calling
        proband_sample = proband_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        pat_sample = pat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        mat_sample = mat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}

        // Prepare the reference channel
        ref = Channel.fromPath(params.ref, checkIfExists: true)
        index_ref = index_ref_fai(ref)
        ref_index = index_ref.reference_index
        // Build ref cache for CRAM steps that do not take a reference
        cram_cache(ref)
        ref_cache = cram_cache.out.ref_cache
        ref_path = cram_cache.out.ref_path
        // canonical ref and BAM channels to pass around to all processes
        ref_channel = ref.concat(ref_index).concat(ref_cache).concat(ref_path).buffer(size: 4)

        // Check Env - Clair3 trio version
        // to get contigs and tmp folder
        // Use clair3 nova model for trio part of the workflow
        checkEnv_trio(proband_sample, pat_sample, mat_sample, snp_bed, ref_channel, chromosome_codes, clair3_nova_model)
        contigs = checkEnv_trio.out.contigs_file.splitText() { it.trim() }

        // use clair3 split beds if BED was provided
        if (params.bed) {
            using_user_bed = true
            split_beds = checkEnv_trio.out.split_beds
        }
        else {
            using_user_bed = false
            split_beds = Channel.from(OPTIONAL_FILE).collect()
        }


        // run clair3 on each sample
        // wf-human-snp also been updated to match clair3 version in nova workflow - may decide to swap back to snp but this seems to lead to different results
        // Use normal clair3 model for pile up part of the workflow (which is largly just the clair3 pileup part)
        proband_snp = snp_proband(proband_sample, snp_bed, ref_channel, clair3_model, extensions, using_user_bed, chromosome_codes)
        pat_snp = snp_pat(pat_sample, snp_bed, ref_channel, clair3_model, extensions, using_user_bed, chromosome_codes)
        mat_snp = snp_mat(mat_sample, snp_bed, ref_channel, clair3_model, extensions, using_user_bed, chromosome_codes)

        // Select candidates per contig
        selectCandidates_trio(contigs, proband_snp.vcf_pileup, pat_snp.vcf_pileup, mat_snp.vcf_pileup)
    
        // Gather all phased bams and combine with the candidate bed files.
        gather_bams = proband_snp.haplotagged_bams
            .join(pat_snp.haplotagged_bams)
            .join(mat_snp.haplotagged_bams)
            .combine(selectCandidates_trio.out.candidate_bed.transpose(), by:0)
        // Call var bam trio per chunk in candidate bed file
        callVarBam_trio(gather_bams.combine(ref_channel).combine(checkEnv_trio.out.tmp).combine(clair3_nova_model))

        // Mix the calls output to create one channel
        calls = callVarBam_trio.out.proband_calls.mix(callVarBam_trio.out.pat_calls, callVarBam_trio.out.mat_calls)
   
        // Group all contig.vcfs per sample so they can be saved in a per-sample directory as required for sortVcf
        all_contig_calls = calls.map{ contig, meta, vcf -> [meta, vcf]}.groupTuple()
        sortVcf_Trio(all_contig_calls.combine(ref_channel).combine(checkEnv_trio.out.tmp))

        // Combine pileup and trio vcfs for merging
        vcf_sets = sortVcf_Trio.out.trio_pileup.combine(proband_snp.vcf_pileup.mix(pat_snp.vcf_pileup, mat_snp.vcf_pileup), by: 0)

        // GVCFs
        aggregated_gvcf = proband_snp.agg_gvcf.mix(pat_snp.agg_gvcf, mat_snp.agg_gvcf)

        // Merge vcfs per contig
        contigs 
        | combine(
            vcf_sets
           | combine(ref_channel)
           | join(aggregated_gvcf)
        )
        | combine(selectCandidates_trio.out.candidate_bed_folder, by: 0)
        | mergeVcf_Trio

        // Group merged vcfs by sample for sorting. Contig is needed for phase_trio below.
        grouped_vcf = mergeVcf_Trio.out.merged_vcf
            .map{ meta, contig, vcf -> [meta, vcf]}
            .groupTuple()

        sortmergedVcf_Trio(
            grouped_vcf
            .combine(ref_channel)
            .combine(checkEnv_trio.out.tmp))
         
        mergeVcf_Trio.out.merged_gvcf
            | map{ meta, contig, vcf -> [meta, vcf]}
            | groupTuple()
            | combine(ref_channel)
            | combine(checkEnv_trio.out.tmp)
            | sortmergedGVcf_Trio

        // Create tuples to name bams, for combining with matching vcfs for phasing
        pat_bam_named = pat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        mat_bam_named = mat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        proband_bam_named = proband_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}

        // Phasing
        // TODO parameterise
        if (params.phase_trio){
            bams_to_phase = pat_bam_named.mix(mat_bam_named, proband_bam_named)
            mergeVcf_Trio.out.merged_vcf | combine(bams_to_phase, by:0) | combine(ref_channel) | phase_trio
            grouped_phased_trio = phase_trio.out.phased.map{ meta, contig, vcf -> [meta, vcf]}.groupTuple()
            grouped_phased_trio | combine(ref_channel) | combine(checkEnv_trio.out.tmp) | sortphased_Trio
            phase_trio.out.phased | combine(bams_to_phase, by:0) | combine(ref_channel) | haplotag_trio
            cat_haplotagged_contigs(haplotag_trio.out.haplotag_bam.groupTuple().combine(ref_channel), extensions)
            hap_bam = cat_haplotagged_contigs.out.merged_xam.flatten()
            phased_vcfs = sortphased_Trio.out.final_vcf.flatMap{ meta, vcf, tbi -> [vcf, tbi]}
        } else {
            hap_bam = Channel.empty()
            phased_vcfs = Channel.empty()
        }

        // Prepare outputs
        vcfs = sortmergedVcf_Trio.out.sorted_trio.flatMap{ meta, vcf, tbi -> [vcf, tbi] }
        gvcfs = sortmergedGVcf_Trio.out.sorted_trio.flatMap{ meta, vcf, tbi -> [vcf, tbi] }


        // Prepare per contig VCFS for GLNEXUS
        //sort gvcfs per contig
        sorted_per_contig = sortgvcf_Trio_contig(mergeVcf_Trio.out.merged_gvcf.combine(ref_channel)
            .combine(checkEnv_trio.out.tmp))
        contig_vcfs_branched = sorted_per_contig
        | branch{
            meta, contig, gvcf, tbi ->
            proband: meta.alias == "${params.proband_sample_name}"
            paternal: meta.alias == "${params.pat_sample_name}"
            maternal: meta.alias == "${params.mat_sample_name}"
        }
        proband_contig_gvcf = contig_vcfs_branched.proband.map{ meta, contig, gvcf, tbi-> [contig, gvcf]}
        paternal_contig_gvcf = contig_vcfs_branched.paternal.map{ meta, contig, gvcf, tbi-> [contig, gvcf]}
        maternal_contig_gvcf = contig_vcfs_branched.maternal.map{ meta, contig, gvcf, tbi-> [contig, gvcf]}
        grouped_gvcf_per_contig = proband_contig_gvcf
        | join(paternal_contig_gvcf)
        | join(maternal_contig_gvcf)

        proband_contig_vcf = mergeVcf_Trio.out.merged_vcf.filter { meta, vcf, tbi ->
            meta.alias == "${params.proband_sample_name}"}.map{ meta, contig, vcf -> [contig, vcf]}
        
        // Run Glnexus for joint variant calling - considering all of the trio simultaniously 
        gln = glnexus(
            grouped_gvcf_per_contig
            | join(proband_contig_vcf)
            | combine(fam_id)
            | combine(gl_conf)
        )
        merged_sorted_vcf = gln | groupTuple | concat_vcfs
        // RTG tools on final vcf to checks a multi-sample VCF file for variant calls which do not follow Mendelian
        //inheritance, and compute aggregate sample concordance
        rtg = rtgTools(merged_sorted_vcf.final_vcf, ref, ped_file)

        // Prepare outputs
        rtg_summary = rtg.summary.map{ family_name, sum -> sum }
        gl_vcf = merged_sorted_vcf.final_vcf.flatMap{ family_name, vcf, tbi -> [vcf, tbi]}

    emit:
        vcf = vcfs
        gvcf = gvcfs
        phased_vcf = phased_vcfs
        haplotagged_bam = hap_bam
        phased_vcf = phased_vcfs
        haplotagged_bam = hap_bam
        multi_vcf = gl_vcf
        rtg_sum = rtg_summary

}