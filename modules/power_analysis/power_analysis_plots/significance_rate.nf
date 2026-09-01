process PLOT_SIGNIFICANCE_RATE {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'

    // Take the input:
    // - stats_file (.csv): path to the summary statistics table used to compute the significance rate
    // - alpha_val (float): significance threshold used to flag a test as significant
    input:
    path stats_file
    val  alpha_val

    // Specify the output of the process
    output:
    path "*.png"

    script:
    """
    plot_significance_rate.py \\
        --csv_file ${stats_file} \\
        --alpha ${alpha_val} \\
        --output_plot "significance_rate_p${alpha_val}.png"
    """
}
