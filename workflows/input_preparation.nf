#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/
include { VALIDATE_INPUTS; BUILD_SAMPLE_REGISTRY } from '../modules/input_management/mutational_matrix/input_validation.nf'
include { GET_SAMPLE_COL; MERGE_SAMPLES } from '../modules/input_management/mutational_matrix/merge_samples.nf'

include { RECONSTRUCT_MUTATIONAL_MATRIX } from '../modules/input_management/reconstruction/reconstruct_mutational_matrix.nf'

include { ALIGN_MUTATION_TYPES } from '../modules/input_management/align_mutation_type_col.nf'

/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/

workflow IP_MM {

    take:
    metadata_path
    metadata_delim
    samplesheet_path
    samplesheet_delim

    main:

    VALIDATE_INPUTS(
        metadata_path,
        metadata_delim,
        samplesheet_path,
        samplesheet_delim
    )

    ch_metadata    = VALIDATE_INPUTS.out.metadata
    ch_samplesheet = VALIDATE_INPUTS.out.samplesheet

    BUILD_SAMPLE_REGISTRY(
        ch_metadata,
        ch_samplesheet
    )

    ch_registry = BUILD_SAMPLE_REGISTRY.out.registry

    GET_SAMPLE_COL(
        ch_registry
    )

    ALIGN_MUTATION_TYPES(
        GET_SAMPLE_COL.out.sample_matrix
            .map { meta, matrix_file -> matrix_file }
    )

    ch_all_matrices = ALIGN_MUTATION_TYPES.out.aligned_matrix
        .collect()

    MERGE_SAMPLES(
        ch_all_matrices
    )

    emit:
    merged_matrix = MERGE_SAMPLES.out.merged_matrix

}

workflow IP_TR {

    take:
    metadata_path
    metadata_delim
    samplesheet_path
    samplesheet_delim

    main:

    // Validate inputs and generate channels
    VALIDATE_INPUTS(
        metadata_path,
        metadata_delim,
        samplesheet_path,
        samplesheet_delim
    )

    // Access the emitted channels
    ch_reconstruction_files = VALIDATE_INPUTS.out.samplesheet

    ch_activities = ch_reconstruction_files.map { files ->
        tuple(
            files[0][0],
            files[0][1]
        )
    }

    ch_signatures = ch_reconstruction_files.map { files ->
        tuple(
            files[1][0],
            files[1][1]
        )
    }
        
    RECONSTRUCT_MUTATIONAL_MATRIX(
        ch_signatures,
        ch_activities
    )

    ch_reconstructed_matrix = RECONSTRUCT_MUTATIONAL_MATRIX.out.reconstructed_matrix

    ALIGN_MUTATION_TYPES(
        ch_reconstructed_matrix
    )

    emit:
    mutational_matrix = ALIGN_MUTATION_TYPES.out.aligned_matrix.collect()

}
