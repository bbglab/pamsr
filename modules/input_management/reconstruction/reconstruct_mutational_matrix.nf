process RECONSTRUCT_MUTATIONAL_MATRIX {

    // Specify the process tag
    tag "reconstruction"

    // Specify the cpus that can be used by the process
    cpus 1

    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Publish the results
    publishDir "${params.outdir}/${params.project_name}/reconstruction", mode: 'copy'

    // Take the input:
    // - signatures (.tsv/.csv): path to the signatures matrix (MutationType x Signature) plus its delimiter ("csv"/"tsv")
    // - activities (.tsv/.csv): path to the activities/exposures matrix (Sample x Signature) plus its delimiter ("csv"/"tsv")
    input:
    tuple path(signatures), val(signatures_delim)
    tuple path(activities), val(activities_delim)

    // Specify the output of the process and emit it
    output:
    path "reconstructed_mutational_matrix.tsv", emit: reconstructed_matrix

    script:

    // Resolve the separator character to use when reading each input file, based on its declared delimiter
    def signatures_sep = signatures_delim == "csv" ? "," : "\t"
    def activities_sep = activities_delim == "csv" ? "," : "\t"

    """
    python3 - <<'PY'

    import pandas as pd
    import numpy as np

    # ============================================================
    # Delimiters
    # ============================================================

    signatures_sep = "${signatures_sep}"
    activities_sep = "${activities_sep}"

    # ============================================================
    # Load matrices
    # ============================================================

    df_prob = pd.read_csv(
        "${signatures}",
        sep=signatures_sep
    )

    df_act = pd.read_csv(
        "${activities}",
        sep=activities_sep
    )

    # ============================================================
    # P: MutationType x Signature
    # ============================================================

    P = df_prob.set_index("MutationType")

    P.columns = P.columns.str.strip()

    # ============================================================
    # E: Sample x Signature
    # ============================================================

    E = df_act.set_index("Samples")

    E.columns = E.columns.str.strip()

    # ============================================================
    # Eliminate target signature
    # ============================================================

    target_signature = "${params.target_signature_eliminate_reconstruction}"

    if target_signature:
        print(f"Eliminating target signature: {target_signature}")

        if target_signature in P.columns:
            P = P.drop(columns=[target_signature])
            print(f"  Removed {target_signature} from P")

        if target_signature in E.columns:
            E = E.drop(columns=[target_signature])
            print(f"  Removed {target_signature} from E")
            
    # ============================================================
    # Check duplicate signature names
    # ============================================================

    dupP = P.columns[P.columns.duplicated()].unique().tolist()
    dupE = E.columns[E.columns.duplicated()].unique().tolist()

    if dupP:
        raise ValueError(
            f"Duplicate signature columns in P: {dupP}"
        )

    if dupE:
        raise ValueError(
            f"Duplicate signature columns in E: {dupE}"
        )

    # ============================================================
    # Align signatures by name
    # ============================================================

    common_sigs = P.columns.intersection(E.columns)

    only_in_P = P.columns.difference(E.columns)
    only_in_E = E.columns.difference(P.columns)

    if len(only_in_P):
        print(
            "Signatures only in P and therefore dropped:",
            list(only_in_P)
        )

    if len(only_in_E):
        print(
            "Signatures only in E and therefore dropped:",
            list(only_in_E)
        )

    if len(common_sigs) == 0:
        raise ValueError(
            "No common signatures found between P and E."
        )

    P_aligned = P[common_sigs]
    E_aligned = E[common_sigs]

    # ============================================================
    # Convert to numeric
    # ============================================================

    P_aligned = P_aligned.apply(
        pd.to_numeric,
        errors="raise"
    )

    E_aligned = E_aligned.apply(
        pd.to_numeric,
        errors="raise"
    )

    # ============================================================
    # Matrix multiplication
    #
    # P = MutationType x Signature
    # E = Sample x Signature
    #
    # M = P x E^T
    # ============================================================

    M = (
        P_aligned.to_numpy()
        @ E_aligned.to_numpy().T
    )

    # Reconstructed mutation counts
    M = np.rint(M).astype(int)

    # ============================================================
    # Create output DataFrame
    # ============================================================

    M_df = pd.DataFrame(
        M,
        index=P_aligned.index,
        columns=E_aligned.index
    )

    M_df.index.name = "MutationType"

    # ============================================================
    # Save
    # ============================================================

    M_df.to_csv(
        "reconstructed_mutational_matrix.tsv",
        sep="\\t"
    )

    print("Reconstructed matrix:")
    print(f"  Mutation types: {M_df.shape[0]}")
    print(f"  Samples:        {M_df.shape[1]}")
    print(f"  Signatures:     {len(common_sigs)}")

    PY
    """
}
