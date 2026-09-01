process PLOT_LOGLH {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'

    // Take the input:
    // - stats_file (.csv): path to the summary statistics table containing log-likelihood values
    // - stat_metric (str): name of the log-likelihood-related statistic/metric to plot
    input:
    path stats_file
    val  stat_metric

    // Specify the output of the process
    output:
    path "*.png"

    script:
    """
    plot_loglh_comparison.py \\
        --csv_file ${stats_file} \\
        --stat ${stat_metric} \\
        --output_plot "loglh_comparison_${stat_metric}.png"
    """
}
