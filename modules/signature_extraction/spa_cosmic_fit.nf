process SIGPROFILERASSIGNMENT_COSMIC_FIT {
    // Specify the process tag
    tag "${params.project_name}"
    // Specify the cpus that can be used by the process
    cpus params.cpus
    // Load the container
    container 'docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3'
    // Publish the results
    publishDir "${params.outdir}/${params.project_name}", mode: 'copy'
    // Take the input:
    // - matrix (.tsv): path to the matrix with the mutational counts of all samples
    // - reference_signatures (.tsv): path to the file containing the mutational catalog used in the process
    input:
    path(matrix)
    path(reference_signatures)
    // Specify the output of the process and emit it
    output:
    path "sig_profiler_assaignment", emit: spa_output

    script:

    """
    # Create the directory needed for the tool to work
    mkdir -p spa_volume
    
    # Assign the variables needed by the tool
    export SIGPROFILERMATRIXGENERATOR_VOLUME="./spa_volume"
    export SIGPROFILERPLOTTING_VOLUME="./spa_volume"
    export SIGPROFILERASSIGNMENT_VOLUME="./spa_volume"

    # Perform the Signature Refitting
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