process PLOT_SIGNIFICANCE_RATE {
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/signature_detection/plots", mode: 'copy'

    input:
    path stats_csv
    val  alpha_val

    output:
    path "*.png"

    script:
    """
    plot_significance_rate.py \\
        --csv_file ${stats_csv} \\
        --alpha ${alpha_val} \\
        --output_plot "significance_rate_p${alpha_val}.png"
    """
}
