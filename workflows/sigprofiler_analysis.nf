nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: NEXTFLOW MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SIGPROFILERASSIGNMENT_COSMIC_FIT } from '../modules/signature_extraction/spa_cosmic_fit.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Run analysis for the sigprofiler assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow SPA {
    // Take the input:
    // -matrix: mutational matrix (96 channels x n samples)
    take:
    matrix
    // Execute the workflow
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