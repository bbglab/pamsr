process GENERATE_SYNTHETIC_COUNTS {

    // Specify the process tag
    tag "synthetic_counts_duplicate_${duplicate_id}"

    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/${params.project_name}/synthetic_data", mode: 'copy'

    // Take the input:
    // - duplicate_id (int): identifier of the replicate; used to derive a distinct seed
    // - seed (int): seed for this replicate's background draw (make it reproducible but unique)
    // - channel_parameters (.tsv): path to the per-channel Gamma-Poisson parameters (alpha_k, theta_k)
    //   produced by COMPUTE_PARAMS_PER_CHANNEL, used to simulate synthetic mutation counts
    input:
    tuple val(duplicate_id), val(seed)
    path channel_parameters

    // Specify the output of the process and emit it
    output:
    tuple val(duplicate_id), path("synthetic_mutation_matrix_duplicate_${duplicate_id}.tsv"),
        emit: synthetic_matrix

    script:
    """
    python3 - <<'PY'

    import numpy as np
    import pandas as pd
    from scipy.stats import poisson

    # ------------------------------------------------------------
    # Parameters
    # ------------------------------------------------------------

    N = ${params.synthetic_N}
    duplicate_id = ${duplicate_id}
    seed = ${seed}

    # Independent RNG stream per duplicate -> independent background draw
    rng = np.random.default_rng(seed)

    # ------------------------------------------------------------
    # Read channel-specific parameters
    # ------------------------------------------------------------

    params_df = pd.read_csv(
        "${channel_parameters}",
        sep="\\t"
    )

    mutation_types = params_df["Mutation Types"]

    lambdas_k = params_df["N_mean_k"].to_numpy(dtype=float)
    num_channels = len(lambdas_k)

    # ------------------------------------------------------------
    # Sample mutation counts
    #
    # Every duplicate gets its own Poisson draw (this call), so
    # background noise is no longer shared/fixed across duplicates.
    # ------------------------------------------------------------

    lambdas_k = np.tile(lambdas_k, (N, 1))
    counts_k = poisson.rvs(
        mu=lambdas_k,
        random_state=rng
    )

    # ------------------------------------------------------------
    # Convert to mutational matrix format
    # ------------------------------------------------------------

    counts_df = pd.DataFrame(
        counts_k.T,
        columns=[
            f"synthetic_{i + 1}"
            for i in range(N)
        ]
    )

    counts_df.insert(
        0,
        "Mutation Types",
        mutation_types
    )

    # ------------------------------------------------------------
    # Save
    # ------------------------------------------------------------

    output_file = f"synthetic_mutation_matrix_duplicate_{duplicate_id}.tsv"

    counts_df.to_csv(
        output_file,
        sep="\\t",
        index=False
    )

    print(f"Duplicate {duplicate_id} | seed {seed} | Generated synthetic mutation matrix")
    print("Number of channels:", num_channels)
    print("Number of samples:", N)
    print("Shape:", counts_df.shape)

    PY
    """
}