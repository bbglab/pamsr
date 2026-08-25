#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER DEFAULTS
========================================================================================
*/
params.metadata          = "${projectDir}/input/metadata/metadata_signature_extraction.csv"
params.metadata_delim    = "csv"

params.samplesheet       = "${projectDir}/input/samplesheet/samplesheet_signature_extraction.csv"
params.samplesheet_delim = "csv"

params.outdir            = "${projectDir}/output"


/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { INPUT_PREPARATION } from './modules/signature_extraction/input_preparation.nf'


/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/
workflow {

    // 1. Create channels for metadata and samplesheet files with existence validation
    ch_metadata    = Channel.fromPath(params.metadata, checkIfExists: true)
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

    // 2. Call the module process
    INPUT_PREPARATION(
        ch_metadata,
        ch_samplesheet,
        params.metadata_delim,
        params.samplesheet_delim
    )
}