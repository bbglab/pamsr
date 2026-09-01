nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: NEXTFLOW MODULES, SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GSD } from '../subworkflows/generate_synthetic_datasets.nf'
include { SIGNATURE_PRESENCE_TEST } from '../modules/power_analysis/signature_detection.nf'
include { SIGNATURE_PRESENCE_SUMMARY } from '../modules/power_analysis/signature_detection.nf'
include { SIGNATURE_PRESENCE_SUMMARY_MERGE } from '../modules/power_analysis/signature_detection.nf'
include { PLOT_PA } from '../subworkflows/plots_power_analysis.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Run analysis for the power analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PA {
    // Take the input:
    // -matrix: mutational matrix (96 channels x n samples)
    take:
    matrix
    // Execute the workflow
    main:
    // ==========================================
    // Generate the synthetic samples
    // ==========================================
    GSD(
        matrix
    )
    // Extract the synthetic samples
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
    // Clean the input so the following process can iterate propperly
    ch_catalog =
        channel.value(
            file(params.reference_signatures_pa)
        )
    // =================================================
    // Test the presence of the signatures using mSigAct
    // =================================================
    SIGNATURE_PRESENCE_TEST(
        ch_injected_individual,
        ch_catalog
    )
    // Extract the results from mSigAct
    ch_signature_detections = SIGNATURE_PRESENCE_TEST.out.results
    // =================================================
    // Summary the results from mSigAct
    // =================================================   
    SIGNATURE_PRESENCE_SUMMARY(
        ch_signature_detections
    )

    ch_summary_individual = SIGNATURE_PRESENCE_SUMMARY.out.summary
        .map { duplicate_id, iteration, n_mutations, file -> file }
        .collect()

    ch_exposures_with_individual = SIGNATURE_PRESENCE_SUMMARY.out.exposures_with
        .map { duplicate_id, iteration, n_mutations, file -> file }
        .collect()

    ch_exposures_without_individual = SIGNATURE_PRESENCE_SUMMARY.out.exposures_without
        .map { duplicate_id, iteration, n_mutations, file -> file }
        .collect()
    // =================================================
    // Merge the summaries of the results from mSigAct
    // =================================================
    SIGNATURE_PRESENCE_SUMMARY_MERGE(
    ch_summary_individual,
    ch_exposures_with_individual,
    ch_exposures_without_individual
    )

    ch_summary_merged = SIGNATURE_PRESENCE_SUMMARY_MERGE.out.summary

    ch_exposures_with_merged = SIGNATURE_PRESENCE_SUMMARY_MERGE.out.exposures_with

    ch_exposures_without_merged = SIGNATURE_PRESENCE_SUMMARY_MERGE.out.exposures_without
    // =================================================
    // Plot the reslts from mSigAct
    // =================================================
    PLOT_PA(
        ch_summary_merged,
        ch_exposures_with_merged,
        ch_exposures_without_merged
    )

}
