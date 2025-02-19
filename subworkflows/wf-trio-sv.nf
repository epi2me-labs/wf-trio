include {
        sniffles2_joint;
        sniffles2;
        sortVCF;
        filterCalls;
        filterCalls as filterTrioCalls;
} from '../modules/local/wf-trio-sv'

include {
    rtgTools;
    annotate_low_complexity as annotate_low_complexity_individual;
    annotate_low_complexity as annotate_low_complexity_joint;
} from '../modules/local/common'
 
// workflow module
workflow sv_trio {
    take:
        sv_bams
        snp_bed
        ref_channel
        ped_file
        optional_file
        chromosome_codes
    main:
        // SV calling
        bams = sv_bams
        | map {
            meta, xam, xai -> [xam, xai, meta] 
        }
        suffix = "wf_trio_sv"
        
        // tandem repeat BED
        if(params.tr_bed) {
            tr_bed = Channel.fromPath(params.tr_bed, checkIfExists: true)
        } else {
            tr_bed = Channel.fromPath(optional_file, checkIfExists: true)
        }
        // Individual SV calling
        snfs = sniffles2(bams.combine(tr_bed).combine(ref_channel), suffix)
        // Joint SV calling
        snf_trio = sniffles2_joint(snfs.snf.map{ meta, snf -> [meta.family_id, snf]}.groupTuple(), ref_channel, suffix)
        sorted_vcf = sortVCF(snf_trio.vcf, suffix)
        // Filter step using BED if provided and contigs 
        filteredCalls = filterTrioCalls(sorted_vcf.vcf_gz.join(sortVCF.out.vcf_tbi).combine(snp_bed), chromosome_codes, suffix)
        final_sv_vcf = filteredCalls
        | map{ xam_meta, vcf, tbi -> [vcf, tbi]}
        | flatten

        // Also filter individual vcfs using BED
        filtered_vcf = filterCalls(snfs.compressed.combine(snp_bed), chromosome_codes, suffix)
        // Checks multi-sample VCF file for variant calls which do not follow Mendelian inheritance
        rtg = rtgTools(filteredCalls, ref_channel, ped_file, suffix)
        rtg_summary_txt = rtg.summary.map{ family_name, sum -> sum }

        // SV annotation
        // Annotate tandem repeat in joint and individual VCFs for downstream processing.  
        prefixed_individual = filtered_vcf.map {xam_meta, vcf, tbi -> [xam_meta.alias, xam_meta, vcf, tbi, "sv"]}
        prefixed_joint = filteredCalls.map{ xam_meta, vcf, tbi -> [xam_meta.family_id, xam_meta, vcf, tbi, "sv"]}
        annotated_individual = annotate_low_complexity_individual(prefixed_individual, "vcf")
        annotated_joint = annotate_low_complexity_joint(prefixed_joint, "vcf")

emit:
    rtg_summary = rtg_summary_txt
}