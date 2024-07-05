// Note this is a reduced version of wf-human-variation.nf (v2.1.0) as a lot of the processes are not needed here
include {
    make_chunks;
    pileup_variants;
    aggregate_pileup_variants;
    select_het_snps;
    phase_contig;
    post_clair_contig_haplotag;
    aggregate_pileup_gvcf;
} from "../modules/local/wf-human-snp.nf"

include { 
    haploblocks as haploblocks_snp
} from '../modules/local/common.nf'

OPTIONAL_FILE = "$projectDir/data/OPTIONAL_FILE"

// workflow module
workflow snp {
    take:
        bam
        bed
        ref
        model
        extensions
        using_user_bed
        chromosome_codes
    main:

        // Run preliminaries to find contigs and generate regions to process in
        // parallel.
        // > Step 0
        make_chunks(bam, ref, bed, model, chromosome_codes)
        meta = bam.map{ meta, bam, bai -> meta }
        chunks = make_chunks.out.chunks_file
            .splitText(){ 
                cols = (it =~ /(.+)\s(.+)\s(.+)/)[0]
                ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]}
        contigs = make_chunks.out.contigs_file.splitText() { it.trim() }
        cmd_file = make_chunks.out.cmd_file
        // use clair3 split beds if BED was provided
        if (using_user_bed) {
            split_beds = make_chunks.out.split_beds
        }
        else {
            split_beds = Channel.from(OPTIONAL_FILE).collect()
        }
        // Run the "pileup" caller on all chunks and collate results
        // > Step 1 
        pileup_variants(chunks, bam, ref, model, bed, cmd_file, split_beds)
        aggregate_pileup_variants(
            ref, pileup_variants.out.pileup_vcf_chunks.collect(),
            make_chunks.out.contigs_file, cmd_file, meta)

        // Filter collated results to produce per-contig SNPs for phasing.
        // > Step 2
        select_het_snps(
            contigs,
            aggregate_pileup_variants.out.pileup_vcf,
            aggregate_pileup_variants.out.phase_qual)
        
        // > Step 3
        // Perform phasing for each contig.
        // `each` does not work with tuples, so we have to make the product ourselves
        phase_inputs = select_het_snps.out.het_snps_vcf
            .combine(bam).combine(ref)
     
        // > Step 4 (haplotagging is now done at the end of the workflow, rather than here)
        phase_contig(phase_inputs)
        phase_contig.out.phased_bam_and_vcf.set { phased_bam_and_vcf }

        // > Step 5
        // merge and sort all files for all chunks for all contigs
        // gvcf is optional, stuff an empty file in, so we have at least one
        // item to flatten/collect and tthis stage can run.
        gvcfs = pileup_variants.out.pileup_gvcf_chunks
            .flatten()
            .ifEmpty(file(OPTIONAL_FILE))
            .collect()
        // aggregate just pile up gvcfs for downstream
        aggregate_pileup_gvcf(ref, gvcfs, make_chunks.out.contigs_file, meta)
        aggregate_gvcf = aggregate_pileup_gvcf.out.agg_gvcf
        // phased requires haplotagged bam to perform appropriate phasing
        // perform internal phasing only if snp+phase is requested, but not sv.
        // Otherwise use final joint phasing only.

        // Only haplotag phased bam needed for trio
        post_clair_contig_haplotag(phased_bam_and_vcf.combine(ref).combine(meta))
        // intermediate ctg BAMs can flow to STR
        haplotagged_ctg_bams = post_clair_contig_haplotag.out.phased_bam
      
        
    emit:
        haplotagged_bams = haplotagged_ctg_bams
        vcf_pileup = aggregate_pileup_variants.out.pileup_vcf
        agg_gvcf = aggregate_gvcf
}


// Reporting workflow
workflow report_snp {
    take:
        vcf_stats
        clinvar_vcf

    main:

        // reporting
        software_versions = getVersions()
        workflow_params = getParams()

        // Create report
        report = makeReport(
            vcf_stats, software_versions.collect(), workflow_params, clinvar_vcf)

    emit:
        report = report
}
