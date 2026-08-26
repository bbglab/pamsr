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
// parameters for sigprofiler
params.reference_signatures = "${projectDir}/input/msigact/Assignment_Solution_Signatures.txt"
params.genome_assembly      = "GRCh38"
params.sigprofiler_args     = ""
params.cpus                 = 12
params.context_type         = "96"

params.project_name         = "testing_spa"
params.synthetic_N = 100
/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/
include { VALIDATE_INPUTS; BUILD_SAMPLE_REGISTRY } from './modules/input_management/input_validation.nf'
include { GET_SAMPLE_COL; MERGE_SAMPLES } from './modules/input_management/merge_samples.nf'
include {
    COMPUTE_PARAMS_PER_CHANNEL
} from './modules/generate_synth_data/compute_parameters.nf'
include {
    GENERATE_SYNTHETIC_COUNTS
} from './modules/generate_synth_data/generate_synthetic_counts.nf'
/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/
workflow {

    // Validate inputs and generate channels
    VALIDATE_INPUTS (
        params.metadata,
        params.metadata_delim,
        params.samplesheet,
        params.samplesheet_delim
    )

    // Access the emitted channels
    ch_metadata    = VALIDATE_INPUTS.out.metadata
    ch_samplesheet = VALIDATE_INPUTS.out.samplesheet

    BUILD_SAMPLE_REGISTRY (
        ch_metadata,
        ch_samplesheet
    )
    ch_registry = BUILD_SAMPLE_REGISTRY.out.registry

    GET_SAMPLE_COL (
        ch_registry
    )

    // Collect all generated file paths into a single channel emission
    ch_all_matrices = GET_SAMPLE_COL.out.sample_matrix
        .map { meta, matrix_file -> matrix_file }
        .collect()
    
    // Merge into a single CSV DataFrame
    MERGE_SAMPLES (
        ch_all_matrices
    )

    ch_mutational_matrix=MERGE_SAMPLES.out.merged_matrix

    COMPUTE_PARAMS_PER_CHANNEL(
        MERGE_SAMPLES.out.merged_matrix
    )

    ch_channel_parameters=COMPUTE_PARAMS_PER_CHANNEL.out.channel_params

    GENERATE_SYNTHETIC_COUNTS(
        ch_channel_parameters
    )
    GENERATE_SYNTHETIC_COUNTS.out.synthetic_matrix.view()
}