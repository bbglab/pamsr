process GENERATE_SYNTHETIC_COUNTS {

    // Specify the process tag
    tag "synthetic_counts_N${params.synthetic_N}"

    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    publishDir "${params.outdir}/${params.project_name}/synthetic_data", mode: 'copy'

    // Take the input:
    // - channel_parameters (.tsv): path to the per-channel Gamma-Poisson parameters (alpha_k, theta_k)
    //   produced by COMPUTE_PARAMS_PER_CHANNEL, used to simulate synthetic mutation counts
    input:
    path channel_parameters

    // Specify the output of the process and emit it
    output:
    path "synthetic_mutation_matrix.tsv", emit: synthetic_matrix

    script:
    """
    python3 - <<'PY'

    import numpy as np
    import pandas as pd
    from scipy.stats import gamma, poisson

    # ------------------------------------------------------------
    # Parameters
    # ------------------------------------------------------------

    N = ${params.synthetic_N}

    # Reproducible random number generator
    rng = np.random.default_rng(42)

    # ------------------------------------------------------------
    # Read channel-specific parameters
    # ------------------------------------------------------------

    params_df = pd.read_csv(
        "${channel_parameters}",
        sep="\\t"
    )

    mutation_types = params_df["Mutation Types"]

    alpha_k = params_df["alpha_k"].to_numpy(dtype=float)
    theta_k = params_df["theta_k"].to_numpy(dtype=float)

    num_channels = len(alpha_k)

    # ------------------------------------------------------------
    # Sample latent lambda values
    #
    # Shape:
    #   (N, K)
    #
    # Each row = one synthetic sample
    # Each column = one mutational channel
    # ------------------------------------------------------------

    # lambdas_k = gamma.rvs(
    #     a=alpha_k,
    #     scale=theta_k,
    #     size=(N, num_channels),
    #     random_state=rng
    # )

    # ------------------------------------------------------------
    # Sample mutation counts
    # ------------------------------------------------------------
    lambdas_k=params_df["N_mean_k"].to_numpy(dtype=float)
    lambdas_k = np.tile(lambdas_k, (N, 1))
    counts_k = poisson.rvs(
        mu=lambdas_k,
        random_state=rng
    )

    # ------------------------------------------------------------
    # Convert to mutational matrix format
    #
    # Current shape:
    #   N x K
    #
    # Required shape:
    #   K x N
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

    counts_df.to_csv(
        "synthetic_mutation_matrix.tsv",
        sep="\\t",
        index=False
    )

    print("Generated synthetic mutation matrix")
    print("Number of channels:", num_channels)
    print("Number of samples:", N)
    print("Shape:", counts_df.shape)

    PY
    """
}
