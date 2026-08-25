import numpy as np
from scipy.stats import gamma, poisson

# ==========================================
# Hierarchical Model Estimators & Sampling
# ==========================================

def estimate_alpha_hierarchical(group):
    """Computes the Method of Moments shape parameter (alpha) for the total

    mutation burden under a hierarchical Gamma-Poisson model.

    The shape parameter alpha governs the underlying Gamma distribution of the
    latent sample-wide mutation rate lambda. Under the Method of Moments for a
    Negative Binomial distribution, alpha is estimated as:

        alpha = (mean_N)^2 / (var_N - mean_N)

    Parameters
    ----------
    group : MutationGroup
        Object instance containing computed sample estimators, specifically:
        - group.N_mean : float
            Empirical mean of total mutation burden across samples.
        - group.s_N_squared : float
            Empirical variance of total mutation burden across samples.

    Returns
    -------
    alpha : float
        Estimated Gamma shape parameter (alpha > 0).

    Raises
    ------
    ValueError
        If empirical variance is less than or equal to the empirical mean
        (var_N <= mean_N), indicating an absence of overdispersion required
        for Negative Binomial modeling.
    """
    mean = group.N_mean
    var = group.s_N_squared

    if np.any(var <= mean):
        raise ValueError(
            "Negative Binomial cannot be estimated because "
            "variance must be strictly larger than the mean."
        )

    alpha = (mean**2) / (var - mean)

    return alpha


def estimate_theta_hierarchical(group):
    """Computes the Method of Moments scale parameter (theta) for the total

    mutation burden under a hierarchical Gamma-Poisson model.

    The scale parameter theta regulates the overdispersion scale of the latent
    sample-wide mutation rate lambda. Under the Method of Moments for a
    Negative Binomial distribution, theta is estimated as:

        theta = (var_N / mean_N) - 1

    Parameters
    ----------
    group : MutationGroup
        Object instance containing computed sample estimators, specifically:
        - group.N_mean : float
            Empirical mean of total mutation burden across samples.
        - group.s_N_squared : float
            Empirical variance of total mutation burden across samples.

    Returns
    -------
    theta : float
        Estimated Gamma scale parameter (theta > 0).

    Raises
    ------
    ValueError
        If empirical variance is less than or equal to the empirical mean
        (var_N <= mean_N), indicating an absence of overdispersion required
        for Negative Binomial modeling.
    """
    mean = group.N_mean
    var = group.s_N_squared

    if np.any(var <= mean):
        raise ValueError(
            "Negative Binomial cannot be estimated because "
            "variance must be strictly larger than the mean."
        )

    theta = (var / mean) - 1.0

    return theta

def sample_mutation_burden_per_sample(group, size=1, random_state=None, return_lambda=False):
    """Generates synthetic total mutation burdens using a hierarchical Gamma-

    Poisson construction (Negative Binomial marginal).

    Sampling proceeds in two hierarchical stages:
    1. Sample latent total mutation rates:
       lambda ~ Gamma(shape=alpha, scale=theta)
    2. Sample discrete total mutation counts given latent rates:
       N ~ Poisson(mu=lambda)

    Parameters
    ----------
    group : MutationGroup
        MutationGroup object containing precomputed estimators for alpha and
        theta derivation.
    size : int, default=1
        Number of synthetic sample burdens to generate.
    random_state : int, np.random.Generator, or None, default=None
        Seed or random number generator instance for reproducible sampling.

    Returns
    -------
    counts : np.ndarray or int
        Synthetic mutation burden(s) sampled from the distribution.
        Shape is (size,) if size > 1, or scalar/1D array depending on SciPy's
        poisson output for size=1.
    """
    alpha = estimate_alpha_hierarchical(group)
    theta = estimate_theta_hierarchical(group)

    # Sample mutation rates per generated sample
    lambdas = gamma.rvs(
        a=alpha, scale=theta, size=size, random_state=random_state
    )

    # Sample total mutation counts
    counts = poisson.rvs(mu=lambdas, random_state=random_state)

    if return_lambda==False:
        return counts
    elif return_lambda==True:
        return counts,lambdas

# ==========================================
# Non-Hierarchical Model Estimators & Sampling
# ==========================================

def estimate_alpha_non_hierarchical(group):
    """Computes channel-specific Method of Moments shape parameters (alpha_k)

    for the non-hierarchical Channel-Specific Mutational Burden (CSMB) model.

    Each channel k is modeled as an independent Gamma-Poisson process.
    The shape vector alpha_k is evaluated as:

        alpha_k = (mean_N_k)^2 / max(var_N_k - mean_N_k, eps)

    Parameters
    ----------
    group : MutationGroup
        Object instance containing channel-specific estimators:
        - group.N_mean_k : pd.Series or np.ndarray
            Mean mutation counts for each channel k across samples.
        - group.s_Nk_squared : pd.Series or np.ndarray
            Variance of mutation counts for each channel k across samples.

    Returns
    -------
    alpha_k : pd.Series or np.ndarray
        Array or Series of estimated Gamma shape parameters for all channels.
    """
    mean_k = group.N_mean_k
    var_k = group.s_Nk_squared

    alpha_k = (mean_k**2) / (var_k - mean_k)

    return alpha_k


def estimate_theta_non_hierarchical(group):
    """Computes channel-specific Method of Moments scale parameters (theta_k)

    for the non-hierarchical Channel-Specific Mutational Burden (CSMB) model.

    Each channel k is modeled as an independent Gamma-Poisson process.
    The scale vector theta_k is evaluated as:

        theta_k = max( (var_N_k / max(mean_N_k, eps)) - 1, eps )

    Parameters
    ----------
    group : MutationGroup
        Object instance containing channel-specific estimators:
        - group.N_mean_k : pd.Series or np.ndarray
            Mean mutation counts for each channel k across samples.
        - group.s_Nk_squared : pd.Series or np.ndarray
            Variance of mutation counts for each channel k across samples.
    eps : float, default=1e-6
        Small positive constant to prevent division by zero for unmutated
        channels and enforce strictly positive scale values (theta_k > 0).

    Returns
    -------
    theta_k : pd.Series or np.ndarray
        Array or Series of estimated Gamma scale parameters for all channels.
    """
    mean_k = group.N_mean_k
    var_k = group.s_Nk_squared

    theta_k = (var_k / mean_k) - 1.0

    return theta_k


def sample_mutation_counts_per_channel(group, size=1, random_state=None, return_lambda=False):
    """Generates synthetic channel-specific mutation count vectors using

    independent Gamma-Poisson processes per channel.

    Sampling proceeds in two stages across all K channels:
    1. Sample channel-specific latent rates:
       lambda_k ~ Gamma(shape=alpha_k, scale=theta_k)
    2. Sample channel-specific discrete counts:
       n_k ~ Poisson(mu=lambda_k)

    Parameters
    ----------
    group : MutationGroup
        MutationGroup object containing precomputed channel estimators.
    size : int, default=1
        Number of synthetic sample count vectors to generate.
    random_state : int, np.random.Generator, or None, default=None
        Seed or random number generator instance for reproducible sampling.

    Returns
    -------
    counts_k : np.ndarray
        Synthetic channel counts matrix or vector.
        - If size == 1: 1D array of shape (K,) containing counts for K channels.
        - If size > 1:  2D array of shape (size, K) where rows are synthetic
          samples and columns are mutation channels.
    """
    alpha_k = estimate_alpha_non_hierarchical(group)
    theta_k = estimate_theta_non_hierarchical(group)

    num_channels = len(alpha_k)

    # Output dimensions based on sample size requests
    sample_shape = (size, num_channels) if size > 1 else num_channels

    # Sample latent baseline mutation rates per channel
    lambdas_k = gamma.rvs(
        a=alpha_k, scale=theta_k, size=sample_shape, random_state=random_state
    )

    # Sample discrete mutation counts per channel
    counts_k = poisson.rvs(mu=lambdas_k, random_state=random_state)

    if return_lambda==False:
        return counts_k
    elif return_lambda==True:
        return counts_k,lambdas_k