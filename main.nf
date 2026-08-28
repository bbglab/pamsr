nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORTS: NEXTFLOW MODULES, WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { IP_MM; IP_TR } from "./workflows/input_preparation.nf"
include { SPA } from "./workflows/sigprofiler_analysis.nf"
include { PA } from "./workflows/power_analysis.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 WORKFLOW: Run main analysis pipeline depending on type of input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    // ===========================
    // Input preparation
    // ===========================
    if (params.input_mode=="mutational_matrix"){
        IP_MM(
        params.metadata_ip_mm,
        params.metadata_delim,
        params.samplesheet_ip_mm,
        params.samplesheet_delim
        )

        ch_mutational_matrix=IP_MM.out.merged_matrix
    }
    else if (params.input_mode == 'tuned_reconstruction') {
        IP_TR(
        params.metadata_ip_tr,
        params.metadata_delim,
        params.samplesheet_ip_tr,
        params.samplesheet_delim
        )

        ch_mutational_matrix=IP_TR.out.mutational_matrix
    }    
    else {

        error """
        Unknown input mode: '${params.input_mode}'

        Valid modes are:
          --input_mode mutational_matrix
          --input_mode tuned_reconstruction
        """
    }
    
    // ===========================
    // Analysis
    // ===========================
    if (params.analysis_mode == 'spa') {

        SPA(
            ch_mutational_matrix
        )

    }
    else if (params.analysis_mode == 'pa') {

        PA(
           ch_mutational_matrix 
        )
        // PA.out.pepe.view()
        // PA.out.a.view()
        // PA.out.b.view()
        // PA.out.c.view()

    }
    else {
        error """
        Unknown pipeline analysis_mode: '${params.analysis_mode}'

        Valid modes are:
          --analysis_mode spa
          --analysis_mode pa
        """
    }
}