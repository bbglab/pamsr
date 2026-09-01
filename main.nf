nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { IP_MM; IP_TR } from "./workflows/input_preparation.nf"
include { SPA } from "./workflows/sigprofiler_analysis.nf"
include { PA } from "./workflows/power_analysis.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Run main analysis pipelines depending on the specified input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    // ==========================================
    // Input preparation from Mutational Matrices
    // ==========================================
    if (params.input_mode=="mutational_matrix"){
        IP_MM(
        params.metadata_ip_mm,
        params.metadata_delim,
        params.samplesheet_ip_mm,
        params.samplesheet_delim
        )

        ch_mutational_matrix=IP_MM.out.merged_matrix
    }
    // =============================================
    // Input preparation from a Tuned Reconstruction
    // =============================================
    else if (params.input_mode == 'tuned_reconstruction') {
        IP_TR(
        params.metadata_ip_tr,
        params.metadata_delim,
        params.samplesheet_ip_tr,
        params.samplesheet_delim
        )

        ch_mutational_matrix=IP_TR.out.mutational_matrix
    }
    // Raise an Error if the input mode is unkown
    else {
        error """
        Unknown input mode: '${params.input_mode}'

        Valid modes are:
          --input_mode mutational_matrix
          --input_mode tuned_reconstruction
        """
    }
    
    // ===========================
    // Signature Profiler Analysis
    // ===========================
    if (params.analysis_mode == 'spa') {

        SPA(
            ch_mutational_matrix
        )

    }
    // ===========================
    // Power analysis
    // ===========================
    else if (params.analysis_mode == 'pa') {

        PA(
           ch_mutational_matrix 
        )
    }
    // Raise an Error if the analysis mode is unkown
    else {
        error """
        Unknown pipeline analysis_mode: '${params.analysis_mode}'

        Valid modes are:
          --analysis_mode spa
          --analysis_mode pa
        """
    }
}