include {
        sniffles2_joint;
        sniffles2;
        sortVCF;
        filterCalls;
        filterCalls as filterTrioCalls;
} from '../modules/local/wf-trio-sv'

include {
    rtgTools;
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
        // Also filter individual vcfs
        compressed_vcf = filterCalls(snfs.compressed.combine(snp_bed), chromosome_codes, suffix)
        // Checks multi-sample VCF file for variant calls which do not follow Mendelian inheritance
        rtg = rtgTools(filteredCalls, ref_channel, ped_file, suffix)
        rtg_summary = rtg.summary.map{ family_name, sum -> sum }
emit:
    rtg = rtg_summary
    sv_unmerged = snfs.vcf
    trio_sv_vcf = final_sv_vcf
    compressed_vcf = compressed_vcf
}