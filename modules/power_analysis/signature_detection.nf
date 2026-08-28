process SIGNATURE_PRESENCE_TEST {

    tag "duplicate_${duplicate_id}_step_${iteration}_n_${n_mutations}"

    container "docker.io/gomdomingoa/msigact:v0.1.0"

    publishDir "${params.outdir}/signature_detection/rdata", mode: 'copy'

    input:
    tuple val(duplicate_id), val(iteration), val(n_mutations), path(matrix)
    path catalog

    output:
    tuple val(duplicate_id), val(iteration), val(n_mutations),
          path("${matrix.simpleName}.sig.presence.test.RData"),
          emit: results

    script:

    def output_name = "${matrix.simpleName}.sig.presence.test.RData"

    """
    test_signature_presence.R \
        ${matrix} \
        ${catalog} \
        ${params.target_signature_pa} \
        ${params.cpus} \
        ${output_name}
    """
}

process SIGNATURE_PRESENCE_SUMMARY {

    tag "duplicate_${duplicate_id}_step_${iteration}_n_${n_mutations}"

    container "docker.io/gomdomingoa/msigact:v0.1.0"

    publishDir "${params.outdir}/signature_detection/summary", mode: 'copy'

    input:
    tuple val(duplicate_id), val(iteration), val(n_mutations), path(rdata)

    output:
    tuple val(duplicate_id), val(iteration), val(n_mutations),
          path("${rdata}_summary_statistics.csv"),
          emit: summary

    tuple val(duplicate_id), val(iteration), val(n_mutations),
          path("${rdata}_exposures_with_target.csv"),
          emit: exposures_with

    tuple val(duplicate_id), val(iteration), val(n_mutations),
          path("${rdata}_exposures_without_target.csv"),
          emit: exposures_without

    script:

    def summary_file =
        "${rdata}_summary_statistics.csv"

    def exp_with_file =
        "${rdata}_exposures_with_target.csv"

    def exp_without_file =
        "${rdata}_exposures_without_target.csv"

    """
    summarize_signature_presence.R \
        ${rdata} \
        ${duplicate_id} \
        ${iteration} \
        ${n_mutations} \
        ${summary_file} \
        ${exp_with_file} \
        ${exp_without_file}
    """
}