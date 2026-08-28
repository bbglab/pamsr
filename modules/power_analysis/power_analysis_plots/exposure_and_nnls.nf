process PLOT_EXPOSURE_NNLS {
    container "docker.io/gomdomingoa/gsd:v0.1.0"
    
    publishDir "${params.outdir}/signature_detection/plots", mode: 'copy'
    
    input:
    path cosmic_tsv
    path exposures_csv
    val  target_signature
    val  bg_signatures // e.g. ["SBS1", "SBS5", "SBS40"]

    output:
    path "*.png"

    script:
    """
    plot_exposure_and_nnls.py \\
        --cosmic_file ${cosmic_tsv} \\
        --exposures_file ${exposures_csv} \\
        --target_sig ${target_signature} \\
        --bg_sigs ${bg_signatures.join(' ')} \\
        --output_plot "exposure_${target_signature}.png"
    """
}