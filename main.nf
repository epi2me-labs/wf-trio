#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended practices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concrete, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2
include {
    lookup_clair3_model;
    lookup_clair3_nova_model;
} from './modules/local/wf-human-snp'
include {
    ingress as proband_ingress; ingress as mat_ingress; ingress as pat_ingress;
} from "./lib/_ingress"
include { snp_trio } from "./subworkflows/wf-trio-snp"
include { sv_trio } from "./subworkflows/wf-trio-sv"
include { getParams } from './lib/common'
include { 
    publish as publish_snp;
    publish as publish_sv;
} from './modules/local/common'
include { concat_vcfs as concat_refined_snp;
          getGenome } from './wf-human-variation/modules/local/common'
include { prepare_reference } from './wf-human-variation/lib/reference'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "wf_common"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    """
}


process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*wf-trio-qc-report.html"
    input:
        tuple val(meta), path(stats), path ("versions/*"), path("params.json")
        path(client_fields)
    output:
        path "*wf-trio-*.html"
    script:
        String report_name = "${meta.alias}.wf-trio-qc-report.html"
        String metadata = new JsonBuilder(meta).toPrettyString()
        String client_fields_args = client_fields.name == OPTIONAL_FILE.name ? "" : "--client_fields $client_fields"
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        --stats $stats/bamstats.readstats.tsv.gz \
        $client_fields_args \
        --params params.json \
        --metadata metadata.json
    """
}


// Creates a new directory named after the sample alias and moves the ingress results
// into it.
process collectIngressResultsInDir {
    label "wf_common"
    input:
        // both inputs might be `OPTIONAL_FILE` --> stage in different sub-directories
        // to avoid name collisions
        tuple val(meta),
            path(reads, stageAs: "reads/*"),
            path(index, stageAs: "index/*"),
            path(stats, stageAs: "stats/*")
    output:
        // use sub-dir to avoid name clashes (in the unlikely event of a sample alias
        // being `reads` or `stats`)
        path "out/*"
    script:
    String outdir = "out/${meta["alias"]}"
    String metaJson = new JsonBuilder(meta).toPrettyString()
    String reads = reads.fileName.name == OPTIONAL_FILE.name ? "" : reads
    String index = index.fileName.name == OPTIONAL_FILE.name ? "" : index
    String stats = stats.fileName.name == OPTIONAL_FILE.name ? "" : stats
    """
    mkdir -p $outdir
    echo '$metaJson' > metamap.json
    mv metamap.json $reads $stats $index $outdir
    """
}


// workflow module
workflow pipeline {
    take:
        reads
    main:
        // fastq_ingress doesn't have the index; add one extra null for compatibility.
        // We do not use variable name as assigning variable name with a tuple
        // not matching (e.g. meta, bam, bai, stats <- [meta, bam, stats]) causes
        // the workflow to crash.
        reads = reads
        .map{
            it.size() == 4 ? it : [it[0], it[1], null, it[2]]
        }

        // Resolve and extract stats files.
        reads.multiMap{ meta, path, index, stats ->
            meta: meta
            stats: stats
        }.set { for_report }
        metadata = for_report.meta.collect()

        client_fields = params.client_fields && file(params.client_fields).exists() ? file(params.client_fields) : OPTIONAL_FILE
        software_versions = getVersions()
        workflow_params = getParams()
        report = makeReport(
           reads.map{meta, path, index, stats -> [meta, stats]}.combine(
            software_versions).combine(workflow_params), client_fields
        )

        // replace `null` with path to optional file
        // If all reads are unaligned path and index may be null
        reads | map { 
            meta, path, index, stats ->
            [ meta, path ?: OPTIONAL_FILE, index ?: OPTIONAL_FILE, stats]
        }
        | collectIngressResultsInDir
    emit:
        ingress_results = collectIngressResultsInDir.out
        workflow_params
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    // demo mutateParam
    if (params.containsKey("mutate_fastq")) {
        CWUtil.mutateParam(params, "fastq", params.mutate_fastq)
    }

    ref = Channel.fromPath(params.ref, checkIfExists: true)
    ref_fp = file(params.ref + ".fai")
    ref_index_fp = Channel.of(ref_fp)
    exts = ["bam", "bai"]

    proband = proband_ingress(
            ref,
            ref_index_fp,
            params.proband_bam,
            exts,
            params.proband_sample_name,
            "proband"
    )
    pat = pat_ingress(
        ref,
        ref_index_fp,
        params.pat_bam,
        exts,
        params.pat_sample_name,
        "pat"
    )
    mat = mat_ingress(
        ref,
        ref_index_fp,
        params.mat_bam,
        exts,
        params.mat_sample_name,
        "mat"
    )
    samples = proband.mix(pat, mat)

    check_genome = getGenome(samples.map{ meta, xam, xai, stats -> [xam, xai, meta]})
    check_genome.map{ genome -> if (genome == "hg19"){
        String input_data_err_msg = '''\
            #####################################################################
            # INPUT DATA PROBLEM
            The genome build detected in the BAM is not compatible with this
            workflow. The workflow does not support hg19.
            See the README for a link to the only supported reference genome hg38/GRCh38
            ################################################################################
            '''.stripIndent()
            error input_data_err_msg
        }
    }

    snp_bed = Channel.fromPath(OPTIONAL_FILE)
    if (params.bed){
        snp_bed = Channel.fromPath(params.bed, checkIfExists: true)
    }

    // Determine basecaller_cfg
    lookup_table = Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
    // attempt to pull out basecaller_cfg from metadata
    def observed_pass_bam = 0
    // Check that basecallers across samples are the same
    metamap_basecaller_cfg = samples
        | map { meta, bam, bai, stats ->
            observed_pass_bam++ // keep count of observed BAMs to guard against emitting basecaller_cfg logging later on failed BAM
            meta["basecall_models"]
        }
        | flatten  // squash lists
    // check returned basecaller list cardinality
    metamap_basecaller_cfg
        | unique
        | count
        | map { int n_models ->
            // n_models of 0 may indicate an empty pass_bam_channel, so
            // we keep a count of observed_pass_bam to determine whether
            // we should handle basecaller_cfg errors
            if (n_models == 0 && observed_pass_bam > 0){
                if (params.override_basecaller_cfg) {
                    log.warn "No basecaller models found in the input alignment header, falling back to the model provided with --override_basecaller_cfg: ${params.override_basecaller_cfg}"
                }
                else {
                    String input_data_err_msg = '''\
                    ################################################################################
                    # INPUT DATA PROBLEM
                    Your input alignment does not indicate the basecall model in the header and you
                    did not provide an alternative with --override_basecaller_cfg.

                    wf-trio requires the basecall model in order to automatically select
                    an appropriate SNP calling model.

                    ## Next steps
                    You must re-run the workflow specifying the basecaller model with the
                    --override_basecaller_cfg option.
                    ################################################################################
                    '''.stripIndent()
                    error input_data_err_msg
                }
            }
            else if (n_models > 1){
                String input_data_err_msg = '''\
                ################################################################################
                # INPUT DATA PROBLEM
                Your input data contains reads basecalled with more than one basecaller model.

                Our workflows automatically select appropriate configuration and models for
                downstream tools for a given basecaller model. This cannot be done reliably when
                reads with different basecaller models are mixed in the same data set.

                ## Next steps
                To use this workflow you must separate your input files, making sure all reads
                have been basecalled with the same basecaller model.
                ################################################################################
                '''.stripIndent()
                error input_data_err_msg
            }
        }

    // use params.override_basecaller_cfg as basecaller_cfg if provided, regardless of what was detected
    // we'll have exploded by now if we have no idea what the config is
    if (params.override_basecaller_cfg) {
        metamap_basecaller_cfg.map {
            log.info "Detected basecaller_model: ${it}"
            log.warn "Overriding basecaller_model: ${params.override_basecaller_cfg}"
        }
        basecaller_cfg = Channel.of(params.override_basecaller_cfg)
    }
    else {
        basecaller_cfg = metamap_basecaller_cfg
            | map { log.info "Detected basecaller_model: ${it}"; it }
            | map { log.info "Using basecaller_model: ${it}"; it }
            | first  // unpack from list
    }


    // map basecalling model to clair3 model
    if(params.clair3_model_path) {
        log.warn "Overriding Clair3 model with ${params.clair3_model_path}."
        clair3_model = Channel.fromPath(params.clair3_model_path, type: "dir", checkIfExists: true)
    }
    else {
        clair3_model = lookup_clair3_model(lookup_table, basecaller_cfg).map {
            log.info "Autoselected Clair3 model: ${it[0]}" // use model name for log message
            it[1] // then just return the path to match the interface above
        }
    }
    // map basecalling model to clair3 nova model
    if(params.clair3_nova_model_path) {
        log.warn "Overriding Clair3 Nova model with ${params.clair3_nova_model_path}."
        clair3_nova_model = Channel.fromPath(params.clair3_nova_model_path, type: "dir", checkIfExists: true)
    }
    else {
        clair3_nova_model = lookup_clair3_nova_model(lookup_table, basecaller_cfg).map {
            log.info "Autoselected Clair3 Nova model: ${it[0]}" // use model name for log message
            it[1] // then just return the path to match the interface above
        }
    }
    fam_id = Channel.of(params.family_id)
    if (params.glnexus_config){
        gl_conf = Channel.fromPath(params.glnexus_config, checkIfExists: true)
    }
    else{
        gl_conf = Channel.fromPath("$projectDir/data/glnexus_conf.yml")
        log.info "No glnexus_config provided using default: $projectDir/data/glnexus_conf.yml"
    }
    ped_file = Channel.fromPath(params.pedigree_file, checkIfExists: true)

    // Prepare the reference channel
    reference = prepare_reference([
        "input_ref": params.ref,
        "output_cache": true,
        "output_mmi": false
    ])
    ref = reference.ref
    ref_index = reference.ref_idx
    ref_cache = reference.ref_cache
    // canonical ref to pass around to all processes
    ref_channel = ref
    | concat(ref_index)
    | concat(ref_cache)
    | flatten
    | buffer(size: 4)

    // Programmatically define chromosome codes.
    ArrayList chromosome_codes = []
    ArrayList chromosomes = [1..22] + ["X", "Y", "M", "MT"]
    for (N in chromosomes.flatten()){
        chromosome_codes += ["chr${N}", "${N}"]
    }


    // Run SNP workflow
    Boolean run_snp = params.snp || (params.sv && params.phased)
    if (run_snp) {
        trio_snp_results = snp_trio(
            proband, pat, mat, snp_bed,
            clair3_model, clair3_nova_model,
            fam_id, gl_conf, ped_file, ref_channel, exts, chromosome_codes)
    
        // Output results
        snp_results = trio_snp_results.rtg_summary | map { [it, null] }

        publish_snp(snp_results)
    }
    // Run SV workflow
    if (params.sv) {
        if (run_snp){
            sv_bams = trio_snp_results.haplotagged_bam
        } else {
            sv_bams = samples.map{ meta, bam, bai, stats -> [meta, bam, bai]}
        }
        trio_sv_results = sv_trio(sv_bams, snp_bed, ref_channel, ped_file, OPTIONAL_FILE, chromosome_codes)
        trio_sv_results.rtg_summary | map { [it, null] } | publish_sv
    }

    // Get stats are standard report
    pipeline(samples)

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
