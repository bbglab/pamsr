nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: NEXTFLOW MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { COMPUTE_PARAMS_PER_CHANNEL } from '../modules/generate_synth_data/compute_parameters.nf'
include { GENERATE_SYNTHETIC_COUNTS } from '../modules/generate_synth_data/generate_synthetic_counts.nf'
include { INJECT_SIGNATURES } from '../modules/power_analysis/signature_injection.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Run pipeline for the generation of synthetic datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow GSD {
    // Take the input:
    // -matrix: mutational matrix (96 channels x n samples)
    take:
    matrix

    main:
    // =============================================
    // Compute the parameters needed for the
    // generation of synthetic samples
    // =============================================
    COMPUTE_PARAMS_PER_CHANNEL(
        matrix
    )

    // .first() turns this into a broadcastable "value channel" so the
    // same channel_parameters file can be reused for every duplicate below
    ch_channel_parameters = COMPUTE_PARAMS_PER_CHANNEL.out.channel_params
    // =============================================
    // Build one (duplicate_id, seed) pair per duplicate,
    // BEFORE simulating the background, so each duplicate
    // gets its own independent Poisson draw.
    // =============================================
    ch_duplicate_seeds = channel
        .fromList( (1..params.n_duplicates.toInteger()).collect { d ->
            tuple(d, params.background_seed.toInteger() + d)
        } )

    // =============================================
    // Generate the synthetic counts (once per duplicate,
    // each with an independent random background)
    // =============================================
    GENERATE_SYNTHETIC_COUNTS(
        ch_duplicate_seeds,
        ch_channel_parameters
    )

    ch_injection_input = GENERATE_SYNTHETIC_COUNTS.out.synthetic_matrix
            
    ch_reference_signatures = channel.value(
        file(params.reference_signatures_gsd)
    )
    // =============================================
    // Generate new samples injecting the mutational
    // signatures at different levels of activities,
    // in an iterative manner
    // =============================================
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