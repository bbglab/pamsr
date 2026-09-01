process PLOT_EXPOSURE_NNLS {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"
    
    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'
    
    // Take the input:
    // - cosmic_tsv (.tsv): path to the COSMIC reference signatures matrix
    // - exposures_csv (.csv): path to the sample exposures table (with the target signature)
    // - target_signature (val): name of the signature under study (e.g. "SBS1")
    // - bg_signatures (val, list): background signatures to include alongside the target in the NNLS model, e.g. ["SBS1", "SBS5", "SBS45"]
    input:
    path cosmic_tsv
    path exposures_csv
    val  target_signature
    val  bg_signatures // e.g. ["SBS1", "SBS5", "SBS40"]

    // Specify the output of the process
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
