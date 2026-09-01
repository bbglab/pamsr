process PLOT_LIMIT_OF_DETECTION {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'

    // Take the input:
    // - stats_file (.csv): path to the summary statistics table used to compute detection limits
    // - alpha_threshold (float): significance threshold used to define the limit of detection
    input:
    path stats_file
    val  alpha_threshold

    // Specify the output of the process
    output:
    path "*.png"

    script:
    """
    plot_limit_of_detection.py \\
        --csv_file ${stats_file} \\
        --alpha ${alpha_threshold} \\
        --output_plot "limit_of_detection_p${alpha_threshold}.png"
    """
}
