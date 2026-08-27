#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/
include { SIGPROFILERASSIGNMENT_COSMIC_FIT } from '../modules/signature_extraction/spa_cosmic_fit.nf'

/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/
workflow SPA {

    take:
    matrix

    main:

    ch_reference_signatures = channel.fromPath(
        params.reference_signatures_spa,
        checkIfExists: true
    )

    SIGPROFILERASSIGNMENT_COSMIC_FIT(
        matrix,
        ch_reference_signatures
    )
    
}