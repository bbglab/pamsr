import numpy as np

from scipy.stats import multinomial

def validate_probabilities(p):
    """
    Validate a multinomial probability vector.

    Parameters
    ----------
    p : numpy.ndarray
        Probability vector.

    Raises
    ------
    ValueError
        If probabilities do not sum to 1 or contain negative values.
    """

    if np.any(p < 0):
        raise ValueError(
            "Multinomial probabilities cannot contain negative values."
        )

    if not np.isclose(np.sum(p), 1):
        raise ValueError(
            f"Multinomial probabilities must sum to 1."
            f"Current sum: {np.sum(p)}"
        )


def synthetic_mutation_counts(
        mutation_burden,
        synthetic_frequencies,
        sample_shape=1,
        random_state=None
    ):
    """
    Generate synthetic mutation counts using a Multinomial distribution.

    Parameters
    ----------
    mutation_burden : int or array-like
        Total number of mutations in the synthetic sample.
        Usually generated from the Negative Binomial model.

    synthetic_frequencies : array-like
        Mutation probability vector generated from the Dirichlet model.
        The probabilities must sum to 1.

        Example:
            [
             p1, p2, ..., p96
            ]

    random_state : int, Generator or None
        Random seed for reproducibility.

    Returns
    -------
    numpy.ndarray
        Synthetic mutation counts for each mutation channel.

        Shape:
            (96,)

    Examples
    --------
    >>> counts = synthetic_mutation_counts(
    ...     mutation_burden=5000,
    ...     synthetic_frequencies=theta,
    ...     random_state=7
    ... )

    """

    p = np.asarray(
        synthetic_frequencies
    )


    validate_probabilities(p)


    counts = multinomial.rvs(
        n=int(mutation_burden),
        p=p,
        size=sample_shape,
        random_state=random_state
        )


    return counts