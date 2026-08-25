// nextflow.enable.dsl=2

params.spectra = "${projectDir}/input/msigact/A_in_data.csv"
params.catalog = "${projectDir}/input/msigact/Assignment_Solution_Signatures.txt"
params.outdir = "${projectDir}/output/msigact"
params.target_sig = "SBS1"
params.cpus = 12

process SIGNATURE_PRESENCE_TEST {

    tag "${params.target_sig}"

    cpus params.cpus

    container "docker.io/gomdomingoa/msigact:v0.1.0"

    input:
    path spectra
    path catalog

    output:
    path "sig.presence.test.out.RData"

    script:
    """
    mkdir -p ${params.outdir}

    test_signature_presence.R \
        ${spectra} \
        ${catalog} \
        ${params.target_sig} \
        ${task.cpus}
    """
}

workflow {

    spectra_ch = Channel.fromPath(
        params.spectra,
        checkIfExists: true
    )

    catalog_ch = Channel.fromPath(
        params.catalog,
        checkIfExists: true
    )

    SIGNATURE_PRESENCE_TEST(
        spectra_ch,
        catalog_ch
    )
}