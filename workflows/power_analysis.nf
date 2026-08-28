#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/

include { GSD } from '../subworkflows/generate_synthetic_datasets.nf'
include { SIGNATURE_PRESENCE_TEST } from '../modules/power_analysis/signature_detection.nf'
include { SIGNATURE_PRESENCE_SUMMARY } from '../modules/power_analysis/signature_detection.nf'

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
    
    ch_injected =
        GSD.out.injected_matrices

    // Convert each group of TSVs into individual tuples
    ch_injected_individual =
        ch_injected.flatMap { duplicate_id, mutation_steps, files ->

            files.withIndex(1).collect { file, index ->

                tuple(
                    duplicate_id,
                    index,
                    mutation_steps[index - 1],
                    file
                )
            }
        }

    ch_catalog =
        channel.value(
            file(params.reference_signatures_pa)
        )

    SIGNATURE_PRESENCE_TEST(
        ch_injected_individual,
        ch_catalog
    )

    ch_signature_detections = SIGNATURE_PRESENCE_TEST.out.results
    
    SIGNATURE_PRESENCE_SUMMARY(
        ch_signature_detections
    )
    
    emit:
    a = SIGNATURE_PRESENCE_SUMMARY.out.summary
    b = SIGNATURE_PRESENCE_SUMMARY.out.exposures_with
    c = SIGNATURE_PRESENCE_SUMMARY.out.exposures_with

}
