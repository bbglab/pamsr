process EXPOSURE_NNLS {
    input:
    path cosmic_csv
    path exposures_csv
    val  target_signature
    val  bg_signatures // e.g. ["SBS1", "SBS5", "SBS40"]

    output:
    path "*.png"

    script:
    """
    plot_exposure_nnls.py \\
        --cosmic_file ${cosmic_csv} \\
        --exposures_file ${exposures_csv} \\
        --target_sig ${target_signature} \\
        --bg_sigs ${bg_signatures.join(' ')} \\
        --output_plot "exposure_${target_signature}.png"
    """
}