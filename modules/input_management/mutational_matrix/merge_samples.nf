process GET_SAMPLE_COL {
    tag "${meta.sample_id}"
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Take the input:
    // - meta (val): a metadata map/object holding, among other fields, matrix_path (path to the source matrix),
    //   sample_id (the sample/column to extract) and delim (delimiter of the source matrix, "csv" or "tsv")
    input:
    val meta

    // Specify the output of the process and emit it
    output:
    tuple val(meta), path("${meta.sample_id}_mut_mat.tsv"), emit: sample_matrix

    script:
    // Resolve the separator character used to read the source matrix, based on meta.delim
    def sep = (meta.delim.toLowerCase() == 'tsv') ? '\\t' : ','
    """
    #!/usr/bin/env python3
    import pandas as pd

    # Read original matrix using meta.matrix_path
    df = pd.read_csv("${meta.matrix_path}", sep="${sep}")

    # Determine label column
    label_col = 'Mutation Types' if 'Mutation Types' in df.columns else df.columns[0]

    # Subset label column + sample column using meta.sample_id
    sample_df = df[[label_col, "${meta.sample_id}"]]

    # Export individual matrix
    sample_df.to_csv("${meta.sample_id}_mut_mat.tsv", sep="\\t", index=False)
    """
}

process MERGE_SAMPLES {
    tag "${params.project_name}"
    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/merged_matrix", mode: 'copy'

    // Take the input:
    // - sample_files (.tsv, collection): the per-sample "*_mut_mat.tsv" matrices produced by GET_SAMPLE_COL, to be merged
    input:
    path sample_files

    // Specify the output of the process and emit it
    output:
    path "merged_mutation_matrix.tsv", emit: merged_matrix

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob

    # Find all sample TSVs
    file_list = sorted(glob.glob("*_mut_mat.tsv"))

    if not file_list:
        raise ValueError("No sample matrix TSV files found to merge.")

    # Read the first file to set the base DataFrame with 'Mutation Types'
    base_df = pd.read_csv(file_list[0],sep="\\t")
    label_col = 'Mutation Types' if 'Mutation Types' in base_df.columns else base_df.columns[0]
    
    # Set the label column as index
    merged_df = base_df.set_index(label_col)

    # Join every subsequent sample file on the index
    for f in file_list[1:]:
        df = pd.read_csv(f,sep="\\t").set_index(label_col)
        merged_df = merged_df.join(df, how='outer')

    # Reset index to restore 'Mutation Types' as the first column
    merged_df.reset_index(inplace=True)

    sample_cols = [c for c in merged_df.columns if c != label_col]

    merged_df[sample_cols] = (
    merged_df[sample_cols]
    .apply(pd.to_numeric, errors="raise")
    .round()
    .astype("int64")
    )

    merged_df.to_csv(
        "merged_mutation_matrix.tsv",
        sep="\\t",
        index=False
    )
    """
}
