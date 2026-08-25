process SIGPROFILERASSIGNMENT_COSMIC_FIT {

    tag "${proj_name}"

    cpus params.cpus

    container "docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3"

    input:
    tuple val(proj_name), path(matrix)
    path(reference_signatures)

    output:
    tuple val(proj_name), path("output_${proj_name}/**/Assignment_Solution_Activities.txt"), emit: sinatures_activities

    script:
    def assembly = task.ext.genome_assembly ?: params.genome_assembly

    """
    mkdir -p spa_volume

    export SIGPROFILERMATRIXGENERATOR_VOLUME="./spa_volume"
    export SIGPROFILERPLOTTING_VOLUME="./spa_volume"
    export SIGPROFILERASSIGNMENT_VOLUME="./spa_volume"


    SigProfilerAssignment cosmic_fit \\
        ${matrix} \\
        output_${proj_name} \\
        --signature_database ${reference_signatures} \\
        --genome_build ${assembly} \\
        --cpu ${task.cpus} \\
        --context_type 96 \\
        --volume spa_volume \\
        ${params.sigprofiler_args}
    """
}