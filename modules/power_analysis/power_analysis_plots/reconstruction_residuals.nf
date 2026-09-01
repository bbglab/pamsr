process PLOT_RECONSTRUCTION_RESIDUALS {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'

    // Take the input:
    // - stats_file (.csv): path to the summary statistics table used to plot reconstruction residuals
    input:
    path stats_file

    // Specify the output of the process
    output:
    path "*.png"

    script:
    """
    plot_reconstruction_residuals.py \\
        --csv_file ${stats_file} \\
        --output_plot "reconstruction_residuals.png"
    """
}
