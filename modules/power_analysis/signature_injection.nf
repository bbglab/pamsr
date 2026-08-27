process INJECT_SIGNATURES {

    tag "duplicate_${duplicate_id}_${params.target_signature_injection}"

    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/injected", mode: 'copy'

    input:
    tuple val(duplicate_id), path(synthetic_matrix)
    path(reference_signatures)
    val(mutation_steps)
    val(target_signature_injection)
    val(seed)

    output:
    path "injected_duplicate_${duplicate_id}_step_*.tsv", emit: injected_matrices

    script:
    """
    python3 - <<'PY'

    import numpy as np
    import pandas as pd

    # ============================================================
    # Parameters
    # ============================================================

    duplicate_id = ${duplicate_id}
    target_signature = "${target_signature_injection}"
    mutation_steps = ${mutation_steps}
    seed = ${seed}

    # ============================================================
    # Input files
    # ============================================================

    df_output = pd.read_csv(
        "${synthetic_matrix}",
        sep="\\t"
    )

    df_signatures = pd.read_csv(
        "${reference_signatures}",
        sep="\\t"
    )

    # ============================================================
    # Mutation labels
    # ============================================================

    mutation_col = df_output.columns[0]

    input_labels = df_output[mutation_col].astype(str)
    signature_labels = df_signatures.iloc[:, 0].astype(str)

    if set(input_labels) != set(signature_labels):
        raise ValueError(
            "Mutation channels do not match between "
            "synthetic matrix and signature matrix."
        )

    # Align signature to synthetic matrix
    df_signatures = (
        df_signatures
        .set_index(df_signatures.columns[0])
        .reindex(input_labels)
    )

    # ============================================================
    # Target signature
    # ============================================================

    if target_signature not in df_signatures.columns:
        raise ValueError(
            f"Target signature '{target_signature}' "
            f"not found in signature matrix."
        )

    p = df_signatures[target_signature].to_numpy(dtype=float)

    if np.any(p < 0):
        raise ValueError(
            f"Signature '{target_signature}' contains negative values."
        )

    if p.sum() <= 0:
        raise ValueError(
            f"Signature '{target_signature}' has zero probability."
        )

    p = p / p.sum()

    # ============================================================
    # Random number generator
    # ============================================================

    rng = np.random.default_rng(
        seed + duplicate_id
    )

    # ============================================================
    # Samples that can receive mutations
    # ============================================================

    sample_columns = list(df_output.columns[1:])

    if len(sample_columns) == 0:
        raise ValueError(
            "Synthetic matrix contains no sample columns."
        )

    # ============================================================
    # Sequential injection
    # ============================================================

    for step, n_mutations in enumerate(
        mutation_steps,
        start=1
    ):

        n_mutations = int(n_mutations)

        # --------------------------------------------------------
        # Each individual mutation is assigned to a random sample
        # --------------------------------------------------------

        for i in range(n_mutations):

            selected_sample = rng.choice(
                sample_columns
            )

            # One mutation drawn from the target signature
            counts = rng.multinomial(
                n=1,
                pvals=p
            )

            df_output.loc[:, selected_sample] += counts

        # --------------------------------------------------------
        # Save the state AFTER this injection step
        # --------------------------------------------------------

        output_file = (
            f"injected_duplicate_{duplicate_id}"
            f"_step_{step}"
            f"_n_{n_mutations}.tsv"
        )

        df_output.to_csv(
            output_file,
            sep="\\t",
            index=False
        )

        print(
            f"Duplicate {duplicate_id} | "
            f"Step {step}/{len(mutation_steps)} | "
            f"Added {n_mutations} mutations | "
            f"Signature {target_signature} | "
            f"Saved {output_file}"
        )

    PY
    """
}