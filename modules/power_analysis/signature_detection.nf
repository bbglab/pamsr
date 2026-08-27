process SIGNATURE_PRESENCE_TEST {

    tag "${matrix.simpleName}"

    container "docker.io/gomdomingoa/msigact:v0.1.0"

    publishDir "${params.outdir}/signature_detection", mode: 'copy'

    input:
    path matrix
    path catalog

    output:
    path "${matrix.simpleName}.sig.presence.test.RData"

    script:
    def output_name = "${matrix.simpleName}.sig.presence.test.RData"

    """
    test_signature_presence.R \\
        ${matrix} \\
        ${catalog} \\
        ${params.target_signature_pa} \\
        ${params.cpus} \\
        ${output_name}
    """
}