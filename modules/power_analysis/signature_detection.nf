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

process SIGNATURE_PRESENCE_SUMMARY_MERGE {

    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/signature_detection/summary", mode: 'copy'

    input:
    path summary_files
    path exposures_with_files
    path exposures_without_files

    output:
    path "merged_summary_statistics.csv", emit: summary
    path "merged_exposures_with_target.csv", emit: exposures_with
    path "merged_exposures_without_target.csv", emit: exposures_without

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob

    # ---------------------------------------------------------
    # Summary statistics
    # ---------------------------------------------------------

    summary_list = sorted(glob.glob("*_summary_statistics.csv"))

    if not summary_list:
        raise ValueError("No summary statistics CSV files found.")

    summary_dfs = []

    for f in summary_list:
        print(f"Reading summary: {f}")
        df = pd.read_csv(f)
        summary_dfs.append(df)

    merged_summary = pd.concat(
        summary_dfs,
        ignore_index=True
    )

    merged_summary.to_csv(
        "merged_summary_statistics.csv",
        index=False
    )


    # ---------------------------------------------------------
    # Exposures with target
    # ---------------------------------------------------------

    with_list = sorted(
        glob.glob("*_exposures_with_target.csv")
    )

    if not with_list:
        raise ValueError("No exposures-with-target CSV files found.")

    with_dfs = []

    for f in with_list:
        print(f"Reading exposures with target: {f}")
        df = pd.read_csv(f)
        with_dfs.append(df)

    merged_with = pd.concat(
        with_dfs,
        ignore_index=True
    )

    merged_with.to_csv(
        "merged_exposures_with_target.csv",
        index=False
    )


    # ---------------------------------------------------------
    # Exposures without target
    # ---------------------------------------------------------

    without_list = sorted(
        glob.glob("*_exposures_without_target.csv")
    )

    if not without_list:
        raise ValueError("No exposures-without-target CSV files found.")

    without_dfs = []

    for f in without_list:
        print(f"Reading exposures without target: {f}")
        without_dfs.append(pd.read_csv(f))

    merged_without = pd.concat(
        without_dfs,
        ignore_index=True
    )

    merged_without.to_csv(
        "merged_exposures_without_target.csv",
        index=False
    )


    print("Merge completed successfully.")
    print(f"Summary files: {len(summary_list)}")
    print(f"With-target files: {len(with_list)}")
    print(f"Without-target files: {len(without_list)}")
    """
}