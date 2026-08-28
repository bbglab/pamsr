process PLOT_LOGLH {
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/signature_detection/plots", mode: 'copy'

    input:
    path stats_csv
    val  stat_metric

    output:
    path "*.png"

    script:
    """
    plot_loglh_comparison.py \\
        --csv_file ${stats_csv} \\
        --stat ${stat_metric} \\
        --output_plot "loglh_comparison_${stat_metric}.png"
    """
}