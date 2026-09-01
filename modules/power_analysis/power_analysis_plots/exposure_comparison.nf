process PLOT_COMPARE_EXPOSURES {
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/signature_detection/plots", mode: 'copy'

    // Take the input:
    // - file_with (.csv): signature exposures CSV computed including the target signature in the presence test model
    // - file_without (.csv): signature exposures CSV computed excluding the target signature in the presence test model
    // - target_signature (str): name of the signature under study (e.g. "SBS1")
    input:
    path file_with
    path file_without
    val  target_signature

    // Specify the output of the process
    output:
    path "*.png"

    script:
    """
    plot_exposure_comparison.py \\
        --file_with ${file_with} \\
        --file_without ${file_without} \\
        --target_sig ${target_signature} \\
        --output_plot "comparison_${target_signature}.png"
    """
}
