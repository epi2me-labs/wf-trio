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
include { lookup_clair3_model } from './modules/local/wf-human-snp'

include {
    ingress as proband_ingress; ingress as mat_ingress; ingress as pat_ingress;
} from "./lib/_ingress"
include {
    trio
} from "./subworkflows/wf-trio"
include {
    getParams;
} from './lib/common'


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
    input:
        val metadata
        path(stats, stageAs: "stats_*")
        path client_fields
        path "versions/*"
        path "params.json"
    output:
        path "wf-trio-*.html"
    script:
        String report_name = "wf-trio-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        String client_fields_args = client_fields.name == OPTIONAL_FILE.name ? "" : "--client_fields $client_fields"
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        --stats $stats \
        $client_fields_args \
        --params params.json \
        --metadata metadata.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wf_common"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
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
           metadata,
           for_report.stats.collect(),
            client_fields,
            software_versions,
            workflow_params
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
        report
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

    snp_bed = Channel.fromPath(OPTIONAL_FILE)
    if (params.bed){
        snp_bed = Channel.fromPath(params.bed, checkIfExists: true)
    }
    if(params.clair3_model_path) {
        log.warn "Overriding Clair3 model with ${params.clair3_model_path}."
        clair3_model = Channel.fromPath(params.clair3_model_path, type: "dir", checkIfExists: true)
    }
    else {
        // map basecalling model to clair3 model
        lookup_table = Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
        clair3_model = lookup_clair3_model(lookup_table, params.basecaller_cfg).map {
            log.info "Autoselected Clair3 model: ${it[0]}" // use model name for log message
            it[1] // then just return the path to match the interface above
        }
    }
    // Run trio pipeline
    trio_results = trio(proband, pat, mat, snp_bed, clair3_model)
    // Get stats are standard report
    samples = proband.mix(pat, mat)
    pipeline(samples)
    // Output results
    pipeline.out.ingress_results
    | map {[ it, null]}
    | concat (
        pipeline.out.report.concat(pipeline.out.workflow_params)
        | map { [it, null] }
    )
    | concat(trio_results.vcf | map { [it, null] })
    | concat(trio_results.gvcf | map { [it, null] })
    | concat(trio_results.phased_vcf | map { [it, null] })
    | concat(trio_results.haplotagged_bam  | map { [it, null] })
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
