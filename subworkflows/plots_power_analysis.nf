#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT SUB-WORKFLOWS & MODULES
========================================================================================
*/

include { PLOT_EXPOSURE_NNLS } from '../modules/power_analysis/power_analysis_plots/exposure_and_nnls.nf'
include { PLOT_COMPARE_EXPOSURES } from '../modules/power_analysis/power_analysis_plots/exposure_comparison.nf'
include { PLOT_LIMIT_OF_DETECTION } from '../modules/power_analysis/power_analysis_plots/limit_of_detection.nf'
include { PLOT_LOGLH } from '../modules/power_analysis/power_analysis_plots/loglh_comparison.nf'
include { PLOT_RECONSTRUCTION_RESIDUALS } from '../modules/power_analysis/power_analysis_plots/reconstruction_residuals.nf'
include { PLOT_SIGNATURE_BLEED } from '../modules/power_analysis/power_analysis_plots/signature_bleed_heatmap.nf'
include { PLOT_SIGNIFICANCE_RATE } from '../modules/power_analysis/power_analysis_plots/significance_rate.nf'

/*
========================================================================================
    WORKFLOW EXECUTION
========================================================================================
*/
workflow PLOT_PA {

    take:
    stats_summary
    exposures_with
    exposures_without

    main:
    PLOT_EXPOSURE_NNLS(
        params.reference_signatures_gsd,
        exposures_with,
        params.target_signature_injection,
        params.background_signatures
    )

    PLOT_COMPARE_EXPOSURES(
        exposures_with,
        exposures_without,
        params.target_signature_compare_exposures
    )

    PLOT_LIMIT_OF_DETECTION(
        stats_summary,
        params.alpha_threshold_log_regression_pa
    )

    PLOT_LOGLH(
        stats_summary,
        params.stat_loglh
    )

    PLOT_RECONSTRUCTION_RESIDUALS(
        stats_summary,
    )

    PLOT_SIGNATURE_BLEED(
        exposures_with,
        exposures_without,
        params.target_signature_pa,
        params.top_n_signatures_heatmap
    )

    PLOT_SIGNIFICANCE_RATE(
        stats_summary,
        params.alpha_threshold_p_value_rate_pa
    )
}