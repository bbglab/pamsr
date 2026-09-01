process ALIGN_MUTATION_TYPES {

    // Specify the process tag
    tag "${matrix.simpleName}"

    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Take the input:
    // - matrix (.tsv): path to a mutational count/probability matrix whose first column holds the mutation type labels
    input:
    path matrix

    // Specify the output of the process and emit it
    output:
    path "aligned_${matrix.simpleName}.tsv", emit: aligned_matrix
    
    script:
    """
    python3 - <<'PY'

    import pandas as pd

    # ------------------------------------------------------------
    # Read input matrix
    # ------------------------------------------------------------

    df = pd.read_csv(
        "${matrix}",
        sep="\\t"
    )

    # ------------------------------------------------------------
    # Standard SBS96 order
    # ------------------------------------------------------------

    major_labels = [
        "C>A",
        "C>G",
        "C>T",
        "T>A",
        "T>C",
        "T>G"
    ]

    flanks = [
        "AA", "AC", "AG", "AT",
        "CA", "CC", "CG", "CT",
        "GA", "GC", "GG", "GT",
        "TA", "TC", "TG", "TT"
    ]

    expected_order = [
        f"{flank[0]}[{subs}]{flank[1]}"
        for subs in major_labels
        for flank in flanks
    ]

    # ------------------------------------------------------------
    # Identify mutation column
    # ------------------------------------------------------------

    mut_col = df.columns[0]

    # ------------------------------------------------------------
    # Align rows to standard SBS96 order
    # ------------------------------------------------------------

    df = (
        df
        .set_index(mut_col)
        .reindex(expected_order)
        .fillna(0.0)
        .reset_index()
    )

    # ------------------------------------------------------------
    # Write aligned matrix
    # ------------------------------------------------------------

    df.to_csv(
        "aligned_${matrix.simpleName}.tsv",
        index=False,
        sep="\\t"
    )

    PY
    """
}
