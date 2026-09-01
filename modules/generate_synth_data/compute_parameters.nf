process COMPUTE_PARAMS_PER_CHANNEL {

    // Specify the process tag
    tag "channel_parameters"

    // Load the container
    container "docker.io/gomdomingoa/gsd:v0.1.0"

    // Take the input:
    // - mutation_matrix (.tsv): path to the merged mutation count matrix (mutation types x samples)
    //   used to estimate per-channel Gamma-Poisson parameters
    input:
    path mutation_matrix

    // Specify the output of the process and emit it
    output:
    path "channel_parameters.tsv", emit: channel_params

    script:
    """
    python3 - <<'PY'

    import pandas as pd

    input_file = "${mutation_matrix}"

    # ------------------------------------------------------------
    # Read merged mutation matrix
    # ------------------------------------------------------------

    df = pd.read_csv(
        input_file,
        sep="\\t"
    )

    # First column = mutation types
    mutation_types = df.iloc[:, 0]

    # Remaining columns = samples
    X = df.iloc[:, 1:].apply(
        pd.to_numeric,
        errors="raise"
    )

    # ------------------------------------------------------------
    # Channel-specific mean and variance
    # ------------------------------------------------------------

    mean_k = X.mean(axis=1)
    var_k = X.var(axis=1, ddof=1)

    # ------------------------------------------------------------
    # Method of Moments estimates
    # ------------------------------------------------------------

    alpha_k = (mean_k ** 2) / (var_k - mean_k)

    theta_k = (var_k / mean_k) - 1.0

    # ------------------------------------------------------------
    # Check parameter constraints
    # ------------------------------------------------------------

    invalid_alpha = alpha_k <= 0
    invalid_theta = theta_k <= 0

    if invalid_alpha.any():

        print("ERROR: The following channels have alpha_k <= 0:")

        print(
            pd.DataFrame({
                "Mutation Types": mutation_types[invalid_alpha],
                "alpha_k": alpha_k[invalid_alpha]
            }).to_string(index=False)
        )

        raise ValueError(
            "All alpha_k parameters must be greater than 0."
        )

    if invalid_theta.any():

        print("ERROR: The following channels have theta_k <= 0:")

        print(
            pd.DataFrame({
                "Mutation Types": mutation_types[invalid_theta],
                "theta_k": theta_k[invalid_theta]
            }).to_string(index=False)
        )

        raise ValueError(
            "All theta_k parameters must be greater than 0."
        )

    # ------------------------------------------------------------
    # Create output
    # ------------------------------------------------------------

    parameters = pd.DataFrame({
        "Mutation Types": mutation_types,
        "N_mean_k": mean_k,
        "s_Nk_squared": var_k,
        "alpha_k": alpha_k,
        "theta_k": theta_k
    })

    parameters.to_csv(
        "channel_parameters.tsv",
        sep="\\t",
        index=False
    )

    print("Channel parameters successfully estimated.")
    print(parameters)

    PY
    """
}
