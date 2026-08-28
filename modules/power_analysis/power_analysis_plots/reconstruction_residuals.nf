process PLOT_RECONSTRUCTION_RESIDUALS {
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/signature_detection/plots", mode: 'copy'

    input:
    path stats_csv

    output:
    path "*.png"

    script:
    """
    plot_reconstruction_residuals.py \\
        --csv_file ${stats_csv} \\
        --output_plot "reconstruction_residuals.png"
    """
}