process PLOT_SIGNATURE_BLEED {
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/signature_detection/plots", mode: 'copy'

    input:
    path file_with_csv
    path file_without_csv
    val  target_sig
    val  top_n_sigs

    output:
    path "*.png"

    script:
    """
    plot_signature_bleed_heatmap.py \\
        --file_with ${file_with_csv} \\
        --file_without ${file_without_csv} \\
        --target_signature ${target_sig} \\
        --top_n ${top_n_sigs} \\
        --output_plot "heatmap_bleed_${target_sig}.png"
    """
}