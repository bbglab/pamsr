#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/

include { GSD } from '../subworkflows/generate_synthetic_datasets.nf'
include { SIGNATURE_PRESENCE_TEST } from '../modules/power_analysis/signature_detection.nf'
/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/
workflow PA {

    take:
    matrix

    main:
    GSD(
        matrix
    )
    
    ch_synthetic_injected_matrices = GSD.out.injected_matrices.flatten()

    ch_catalog =
        channel.value(
            file(params.reference_signatures_pa)
        )

    SIGNATURE_PRESENCE_TEST(
        ch_synthetic_injected_matrices,
        ch_catalog
    )

    emit:
    pepe = ch_synthetic_injected_matrices
}
