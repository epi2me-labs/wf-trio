include { 
    selectCandidates_trio;
    callVarBam_trio;
    sortVcf_Trio;
    mergeVcf_Trio;
    sortmergedVcf_Trio;
    sortmergedGVcf_Trio;
    phase_joint_trio;
    sortphased_Trio as sortJointPhased_Trio;
    haplotag_trio as haplotagJoint_Trio;
    filter_joint;
    cat_haplotagged_contigs as catJoint_haplotagged_contigs;
    checkEnv_trio;
    glnexus;
    sortgvcf_Trio_contig;
    merge_individual_vcfs;
} from '../modules/local/wf-trio-snp'

include {
    concat_vcfs;
    rtgTools;
} from '../modules/local/common'

include {
    snp as snp_proband; snp as snp_mat; snp as snp_pat;
} from "../subworkflows/wf-human-snp"

include {
    concat_vcfs as concat_snp_vcfs;
} from "../wf-human-variation/modules/local/common"

OPTIONAL_FILE = "$projectDir/data/OPTIONAL_FILE"
// workflow module
workflow snp_trio {
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
        ref_channel
        extensions
        chromosome_codes

    main:
        
        

        // Prepare bams required for snp calling
        proband_sample = proband_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        pat_sample = pat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        mat_sample = mat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}

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
         
        mergeVcf_Trio.out.merged_gvcf
            | map{ meta, contig, vcf -> [meta, vcf]}
            | groupTuple()
            | combine(ref_channel)
            | combine(checkEnv_trio.out.tmp)
            | sortmergedGVcf_Trio

        sortmergedVcf_Trio(
                grouped_vcf
                .combine(ref_channel)
                .combine(checkEnv_trio.out.tmp))
                hap_bam = Channel.empty()

        // Prepare outputs
        snp_vcfs = sortmergedVcf_Trio.out.sorted_trio
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
        
        // Run Glnexus for joint variant calling per contig - considering all of the trio simultaniously 
        gln = glnexus(
            grouped_gvcf_per_contig
            | join(proband_contig_vcf)
            | combine(fam_id)
            | combine(gl_conf)
        )
        // Concat per contig joint vcfs
        merged_sorted_vcf = concat_vcfs(gln.map{ family_id, contig, vcf, tbi -> [family_id, vcf]}.groupTuple(), "wf_trio_snp") 
        // RTG tools used on final multi-sample VCF to check for variant calls which do not follow Mendelian
        //inheritance. Compute aggregate sample concordance.
        rtg = rtgTools(merged_sorted_vcf.final_vcf, ref_channel, ped_file, "wf_trio_snp")

        // Prepare outputs
        rtg_summary_txt = rtg.summary.map{ family_name, sum -> sum }
        gl_vcf = merged_sorted_vcf.final_vcf.flatMap{ family_name, vcf, tbi -> [vcf, tbi]}

        // Trio phasing (pedigree joint phasing)
        if (params.phased){
            // Create tuples to name bams, for combining with matching vcfs for phasing
            proband_bam_named = proband_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
            pat_bam_named = pat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
            mat_bam_named = mat_bam.map{ meta, bam, bai, stats -> [meta, bam, bai]}
     
            // Joint pedigree phased
            phase_trio_tuple = glnexus.out.gl_vcf
                    | combine(proband_bam_named)
                    | combine(pat_bam_named)
                    | combine(mat_bam_named)
                    | combine(ref_channel)
                    | combine(ped_file)
            phase_trio_out = phase_joint_trio(phase_trio_tuple)
        
            // Filter joint phased to individual VCFs per contig
            // Order does not matter here 
            bams_to_phase = proband_bam_named
                    | mix(pat_bam_named, mat_bam_named)
            
            // Use phased joint individual VCF's to haplotag bam per contig and haplotagphase
            haplotagJoint_Trio(bams_to_phase
                    | combine(phase_trio_out.joint_phased)
                    | combine(ref_channel))
     
            // Get haplotagphase vcf with phase information added to indels based on haplotagged reads.
            haplotag_phased = haplotagJoint_Trio.out.haplotag_phased_vcf

            // Sort and combine contigs back together
            sortJointPhased_Trio = concat_snp_vcfs(haplotag_phased.groupTuple(), "wf_trio_snp")

            // Merge back to joint VCF
            merge_individual_vcfs(
                params.family_id,
                sortJointPhased_Trio.final_vcf.map{ meta, vcf, tbi -> [vcf, tbi]}.collect())
           
            // Join per contig bam
            catJoint_haplotagged_contigs(
                haplotagJoint_Trio.out.haplotag_bam.map{contig, meta, bam, bai -> [meta, bam]}
                | groupTuple()
                | combine(ref_channel),
                extensions)

            // Get final xam and joint vcf
            trio_hap_bam = catJoint_haplotagged_contigs.out.merged_xam
            vcf_joint = merge_individual_vcfs.out

        } else {
            // If no phasing then just output joint vcf
            vcf_joint = gl_vcf
            trio_hap_bam = Channel.empty()
        }

    emit:
        gvcf = gvcfs
        haplotagged_bam = trio_hap_bam
        rtg_summary = rtg_summary_txt
        contigs = contigs
        joint_vcf = vcf_joint
}