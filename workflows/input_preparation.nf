nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: NEXTFLOW MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { VALIDATE_INPUTS; BUILD_SAMPLE_REGISTRY } from '../modules/input_management/mutational_matrix/input_validation.nf'
include { GET_SAMPLE_COL; MERGE_SAMPLES } from '../modules/input_management/mutational_matrix/merge_samples.nf'

include { RECONSTRUCT_MUTATIONAL_MATRIX } from '../modules/input_management/reconstruction/reconstruct_mutational_matrix.nf'

include { ALIGN_MUTATION_TYPES } from '../modules/input_management/align_mutation_type_col.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Both pipeline for the preparaion of the input, from mutational matrices (MM) or a tuned reconstruction (TR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow IP_MM {
    // Take the input:
    // - metadata_path: path to the metadata file of the samples
    // - metadata_delim: which delimiter is used for the metadata file
    // - samplesheet_path:  path to the sample sheet file that indicates where the input data is located
    // - samplesheet_delim: which delimiter is used for the samplesheet file
    take:
    metadata_path
    metadata_delim
    samplesheet_path
    samplesheet_delim

    main:
    // =============================================
    // Validate the metadata and shample sheet files
    // =============================================
    VALIDATE_INPUTS(
        metadata_path,
        metadata_delim,
        samplesheet_path,
        samplesheet_delim
    )

    ch_metadata    = VALIDATE_INPUTS.out.metadata
    ch_samplesheet = VALIDATE_INPUTS.out.samplesheet
    // ====================================================================
    // Build a registry that identifies in which file is stored each sample
    // ====================================================================
    BUILD_SAMPLE_REGISTRY(
        ch_metadata,
        ch_samplesheet
    )

    ch_registry = BUILD_SAMPLE_REGISTRY.out.registry
    // ======================================
    // Extract each sample as a single column
    // ======================================
    GET_SAMPLE_COL(
        ch_registry
    )
    // ==============================================
    // Align each sample to the desired order of rows
    // ==============================================
    ALIGN_MUTATION_TYPES(
        GET_SAMPLE_COL.out.sample_matrix
            .map { meta, matrix_file -> matrix_file }
    )

    ch_all_matrices = ALIGN_MUTATION_TYPES.out.aligned_matrix
        .collect()
    // ======================================================
    // Merge all the sample columns together in a single file
    // ======================================================
    MERGE_SAMPLES(
        ch_all_matrices
    )

    emit:
    merged_matrix = MERGE_SAMPLES.out.merged_matrix

}

workflow IP_TR {
    // Take the input:
    // - metadata_path: path to the metadata file of the samples
    // - metadata_delim: which delimiter is used for the metadata file
    // - samplesheet_path:  path to the sample sheet file that indicates where the input data is located
    // - samplesheet_delim: which delimiter is used for the samplesheet file
    take:
    metadata_path
    metadata_delim
    samplesheet_path
    samplesheet_delim

    main:
    // =============================================
    // Validate the metadata and shample sheet files
    // =============================================
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
    // =============================================
    // Reconstruct the mutational matrix from the results of a matrix decomposition, eliminating the contribuion of an specific mutational signature
    // =============================================
    RECONSTRUCT_MUTATIONAL_MATRIX(
        ch_signatures,
        ch_activities
    )

    ch_reconstructed_matrix = RECONSTRUCT_MUTATIONAL_MATRIX.out.reconstructed_matrix
    // ==============================================================
    // Align the matrix with all samples to the desired order of rows
    // ==============================================================
    ALIGN_MUTATION_TYPES(
        ch_reconstructed_matrix
    )

    emit:
    mutational_matrix = ALIGN_MUTATION_TYPES.out.aligned_matrix.collect()

}
