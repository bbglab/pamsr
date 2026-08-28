process PLOT_COMPARE_EXPOSURES {
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/signature_detection/plots", mode: 'copy'

    input:
    path file_with_csv
    path file_without_csv
    val  target_signature

    output:
    path "*.png"

    script:
    """
    plot_exposure_comparison.py \\
        --file_with ${file_with_csv} \\
        --file_without ${file_without_csv} \\
        --target_sig ${target_signature} \\
        --output_plot "comparison_${target_signature}.png"
    """
}