/*
 * This sub-workflow acts as a wrapper around `ingress.nf`,
 * extending its functionality to do additional input preprocessing:
 *  - check whether the aligned BAM files match the input reference
 *  - (re-)align the input BAM if unaligned or not matching the reference
 * Additionally for the wf-trio it adds the relationship (proband, maternal, paternal)
 * to meta and reheaders BAM with the relationship
 */
include { xam_ingress } from './ingress.nf'
include { cram_to_bam ;
          minimap2_alignment ;
          check_for_alignment} from "../wf-human-variation/lib/_ingress"
include { reheader_BAM } from '../modules/local/common'


// Create ingress workflow, wrapping functions from `ingress.nf` and
// performing additional processings.
workflow ingress {
    take:
        ref_file
        ref_idx_file
        bam_file_fp
        alignment_exts
        sample_name // TODO: parameterise sample_name in humvar_ingress and re-use
        relationship
    main:
        // load bam as channel
        // We do not want to perform statistics here as we will do them downstream.
        // We also want to keep all unaligned reads.
        ingressed_bam = xam_ingress([
            "input":bam_file_fp,
            "sample":sample_name,
            "sample_sheet":null,
            "analyse_unclassified":true,
            "keep_unaligned": true,
            "stats": true,
            "watch_path": false,
            "per_read_stats": true
        ])
        // Check that we have a single BAM/folder with BAMs in it.
        // by counting how many entries are in the channel.
        // If there are more than 1, then throw an error.
        ingressed_bam
            .count()
            .subscribe { int n_samples -> 
                if (n_samples > 1){
                    error "Too many samples found: (${n_samples}) in ${bam_file_fp}.\nPlease, ensure you provide a single folder with all the BAM files for a single individual."
                }
            }

        // Prepare reference channel
        check_ref = ref_file.combine(ref_idx_file) // don't wait for cram_cache to perform check_for_alignment
        
        // Index the input BAM, and then check if the header matches the reference.
        // Add also if it is a CRAM (for downstream compatibility) and move if they need realignment
        // to meta.
        checked_bam = check_for_alignment(
                check_ref,
                ingressed_bam.map{ meta, xam, xai, stats -> [meta, xam, xai]} 
            ) | 
            map{ to_align, has_maps, meta, xam, xai ->
                [meta + [to_align: to_align != '0', has_mapped_reads: has_maps == '1'], xam, xai]
            }
        // fork BAMs into to_align and noalign subchannels
        checked_bam.branch {
            meta, xam, xai-> 
            to_align: meta.is_unaligned || meta.to_align
            noalign: true
        }.set{alignment_fork}

    
        // call minimap on bams that require (re)alignment
        // then map the result of minimap2_alignment to the canonical (reads, index, meta) tuple
        new_mapped_bams = minimap2_alignment(ref_file, alignment_fork.to_align, alignment_exts).map{
            meta, has_maps, xam, xai -> 
            [meta + [has_mapped_reads: has_maps == '1'], xam, xai]
        }

        bam_channel = alignment_fork.noalign
            .mix(new_mapped_bams)
        
        samples_with_relationship = bam_channel.map{ meta, xam, xai-> [meta + [relationship: relationship, family_id:params.family_id], xam, xai]}   
        reheadered_samples =  reheader_BAM(samples_with_relationship)

        // Add back stats
        samples_with_stats = reheadered_samples.map{ meta, xam, xai -> [meta.alias, meta, xam, xai]}
        .join(ingressed_bam
            .map{ old_meta, xam, xai, stats -> [old_meta.alias, stats]}, by:0)
        .map{ alias, meta, xam, xai, stats -> [meta, xam, xai, stats]}

    emit:
        samples_with_stats
}
