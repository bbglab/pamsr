process SIGPROFILERASSIGNMENT_COSMIC_FIT {

    tag params.project_name

    cpus params.cpus

    container 'docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3'
    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(matrix)
    path(reference_signatures)

    output:
    path "sig_profiler_assaignment", emit: spa_output

    script:

    """
    mkdir -p spa_volume

    export SIGPROFILERMATRIXGENERATOR_VOLUME="./spa_volume"
    export SIGPROFILERPLOTTING_VOLUME="./spa_volume"
    export SIGPROFILERASSIGNMENT_VOLUME="./spa_volume"


    SigProfilerAssignment cosmic_fit \\
        ${matrix} \\
        sig_profiler_assaignment \\
        --signature_database ${reference_signatures} \\
        --genome_build ${params.genome_assembly} \\
        --cpu ${task.cpus} \\
        --context_type 96 \\
        --volume spa_volume \\
        ${params.sigprofiler_args}
    """
}