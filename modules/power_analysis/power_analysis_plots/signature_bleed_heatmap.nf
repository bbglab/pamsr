process PLOT_SIGNATURE_BLEED {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'

    // Take the input:
    // - file_with_target (.csv): signature exposures CSV computed including the target signature in the presence test model
    // - file_without_target (.csv): signature exposures CSV computed excluding the target signature in the presence test model
    // - target_sig (str): name of the signature under study (e.g. "SBS1")
    // - top_n_sigs (int): number of top "bleeding" signatures to show in the heatmap
    input:
    path file_with_target
    path file_without_target
    val  target_sig
    val  top_n_sigs

    // Specify the output of the process
    output:
    path "*.png"

    script:
    """
    plot_signature_bleed_heatmap.py \\
        --file_with ${file_with_target} \\
        --file_without ${file_without_target} \\
        --target_signature ${target_sig} \\
        --top_n ${top_n_sigs} \\
        --output_plot "heatmap_bleed_${target_sig}.png"
    """
}
