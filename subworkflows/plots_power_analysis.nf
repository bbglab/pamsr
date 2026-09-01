nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: NEXTFLOW MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PLOT_EXPOSURE_NNLS } from '../modules/power_analysis/power_analysis_plots/exposure_and_nnls.nf'
include { PLOT_COMPARE_EXPOSURES } from '../modules/power_analysis/power_analysis_plots/exposure_comparison.nf'
include { PLOT_LIMIT_OF_DETECTION } from '../modules/power_analysis/power_analysis_plots/limit_of_detection.nf'
include { PLOT_LOGLH } from '../modules/power_analysis/power_analysis_plots/loglh_comparison.nf'
include { PLOT_RECONSTRUCTION_RESIDUALS } from '../modules/power_analysis/power_analysis_plots/reconstruction_residuals.nf'
include { PLOT_SIGNATURE_BLEED } from '../modules/power_analysis/power_analysis_plots/signature_bleed_heatmap.nf'
include { PLOT_SIGNIFICANCE_RATE } from '../modules/power_analysis/power_analysis_plots/significance_rate.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Run pipeline for plotting the results of the power analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PLOT_PA {
    // Take the input:
    // - stats_summary: file with the statistics of the signature presence test (mSigAct)
    // - exposures_with: signature mutational activities detected when considering the target signature in the model
    // - exposures_without: signature mutational activities detected when NOT considering the target signature in the model 
    take:
    stats_summary
    exposures_with
    exposures_without

    main:
    // ===========================================================
    // Plot the detected exposure of the target signature (x-axis)
    // across the different sample sizes (y-axis) and indicating the
    // NNLS reconstruction error
    // ===========================================================
    PLOT_EXPOSURE_NNLS(
        params.reference_signatures_gsd,
        exposures_with,
        params.target_signature_injection,
        params.background_signatures
    )
    // ===========================================================
    // Plot the detected exposure of the target signature (x-axis)
    // across the different sample sizes (y-axis), comparing the
    // the model that considers the target signature, and the model
    // that does not consider it
    // ===========================================================
    PLOT_COMPARE_EXPOSURES(
        exposures_with,
        exposures_without,
        params.target_signature_compare_exposures
    )
    // ===========================================================
    // Plot of the logistic regression that estimates how much
    // injected signature is needed to precisely detect the 
    // target signature
    // ===========================================================
    PLOT_LIMIT_OF_DETECTION(
        stats_summary,
        params.alpha_threshold_log_regression_pa
    )
    // ===========================================================
    // Plot comparing the log-likelihood of the models (with and
    // without the target signature) for reconstructing the 
    // samples (y-axis) along the different sample sizes (x-axis)
    // ===========================================================
    PLOT_LOGLH(
        stats_summary,
        params.stat_loglh
    )
    // ===========================================================
    // Plot the reconstruction residuals of the total amount of 
    // detected signature activities (y-axis) along the different
    // amounts of signature injections
    // ===========================================================
    PLOT_RECONSTRUCTION_RESIDUALS(
        stats_summary,
    )
    // ===========================================================
    // Plot the heatmap of the different signatures (y-axis) along
    // the different sample sizes (x-axis), coloring each cell
    // based on the signature bleeding
    // ===========================================================
    PLOT_SIGNATURE_BLEED(
        exposures_with,
        exposures_without,
        params.target_signature_pa,
        params.top_n_signatures_heatmap
    )
    // ===========================================================
    // Plot the percentage of samples with a significant detection
    // of the target signature (y-axis) along the different sample
    // sizes (x-axis)
    // ===========================================================
    PLOT_SIGNIFICANCE_RATE(
        stats_summary,
        params.alpha_threshold_p_value_rate_pa
    )
}