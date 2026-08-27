#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/

include { COMPUTE_PARAMS_PER_CHANNEL } from '../modules/generate_synth_data/compute_parameters.nf'
include { GENERATE_SYNTHETIC_COUNTS } from '../modules/generate_synth_data/generate_synthetic_counts.nf'
include { INJECT_SIGNATURES } from '../modules/power_analysis/signature_injection.nf'

/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/
workflow GSD {

    take:
    matrix

    main:
    COMPUTE_PARAMS_PER_CHANNEL(
        matrix
    )

    ch_channel_parameters=COMPUTE_PARAMS_PER_CHANNEL.out.channel_params

    GENERATE_SYNTHETIC_COUNTS(
        ch_channel_parameters
    )

    ch_synthetic = GENERATE_SYNTHETIC_COUNTS.out.synthetic_matrix

    ch_injection_input = ch_synthetic
        .flatMap { synthetic_matrix ->
            (1..params.n_duplicates.toInteger()).collect { duplicate_id ->
                tuple(
                    duplicate_id,
                    synthetic_matrix
                )
            }
        }
            
    ch_reference_signatures = channel.value(
        file(params.reference_signatures_gsd)
    )

    INJECT_SIGNATURES(
        ch_injection_input,
        ch_reference_signatures,
        params.injected_mutations,
        params.target_signature_injection,
        params.injection_seed
    )

    emit:
    injected_matrices = INJECT_SIGNATURES.out.injected_matrices
}