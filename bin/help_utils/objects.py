import numpy as np
import pandas as pd

class MutationGroup:
    """
    Class for all information associated with a single experimental group.

    Parameters
    ----------
    group_id : str
        Identifier of the group.

    samples : pandas.DataFrame
        Mutation count matrix for the samples belonging to this group.
        Rows correspond to mutation channels and columns correspond to samples.

    metadata : pandas.DataFrame
        Metadata corresponding to the samples in the group.
    """

    def __init__(self, group_id, samples, metadata, mutation_types):

        self.group_id = group_id
        self.samples = samples
        self.metadata = metadata
        self.mutation_types = list(mutation_types)

        self.N_g = None
        self.f_kg = None
        self.f_mean_k = None
        self.s_fk_squared = None
        self.N_mean = None
        self.s_N_squared = None

        # Estimators used specifically for the non-hierarchical model (vectors per channel k)
        self.N_mean_k = None
        self.s_Nk_squared = None

    def _validate_non_negative(self, estimator_name, data):
        """Helper to ensure estimators are >= 0."""
        arr = np.asarray(data)
        if np.any(arr < 0):
            raise ValueError(
                f"Validation failed for '{estimator_name}' in group '{self.group_id}': "
                f"Found negative values."
            )

    def __str__(self):
        return (
            f"MutationGroup("
            f"group_id={self.group_id}, "
            f"samples_shape={self.samples.shape}, "
            f"metadata_shape={self.metadata.shape}, "
            f"N_g={self.N_g}, "
            f"f_kg={self.f_kg}, "
            f"f_mean_k={self.f_mean_k}, "
            f"s_fk_squared={self.s_fk_squared}, "
            f"N_mean={self.N_mean}, "
            f"s_N_squared={self.s_N_squared}, "
            f"N_mean_k={self.N_mean_k}, "
            f"s_Nk_squared={self.s_Nk_squared}"
            f")"
        )

    def compute_estimators(self):
        """
        Compute all estimators associated with the mutation group.

        Computes:

        - N_g - Total number of mutations in a given sample g
        - f_kg - Mutation frequency of a given channel in a given sample
        - f_mean_k - Arithmetic mean of the frequency of mutations for a given channel in the whole dataset
        - s_fk_squared - Variance of mutation frequency for a given channel in the whole dataset
        - N_mean - Arithmetic mean of the total number of mutations
        - s_N_squared - Variance of the total number of mutations
        - N_mean_k - Variance of the total number of mutations per channel
        - s_N_squared_k - Variance of the total number of mutations per channel
        """

        X = self.samples
        # Total mutations per sample
        self.N_g = X.sum(axis=0)
        # Mutation frequencies
        self.f_kg = X.div(self.N_g, axis=1)
        # Mean mutation frequency per channel
        self.f_mean_k = self.f_kg.mean(axis=1)
        # Variance of mutation frequency per channel
        self.s_fk_squared = self.f_kg.var(axis=1, ddof=1)
        # Mean of the total mutational burden
        self.N_mean = self.N_g.mean()
        # Variance of the total mutational burden
        self.s_N_squared = self.N_g.var(ddof=1)
        # Mean of the mutational burden per channel
        self.N_mean_k = X.mean(axis=1)
        # Variance of the mutational burden per channel
        self.s_Nk_squared = X.var(axis=1, ddof=1)

        # Validation Step
        estimators = {
            "N_g": self.N_g,
            "f_kg": self.f_kg,
            "f_mean_k": self.f_mean_k,
            "s_fk_squared": self.s_fk_squared,
            "N_mean": self.N_mean,
            "s_N_squared": self.s_N_squared,
            "N_mean_k": self.N_mean_k,
            "s_Nk_squared": self.s_Nk_squared,
        }

        for name, value in estimators.items():
            self._validate_non_negative(name, value)

    def __repr__(self):

        return (
            f"MutationGroup(\n"
            f"  group='{self.group_id}',\n"
            f"  samples={self.samples.shape[1]},\n"
            f"  channels={self.samples.shape[0]},\n"
            f"  N_mean={self.N_mean},\n"
            f"  s_N_squared={self.s_N_squared},\n"
            f"  N_g={'computed' if self.N_g is not None else 'None'},\n"
            f"  f_kg={'computed' if self.f_kg is not None else 'None'},\n"
            f"  f_mean_k={'computed' if self.f_mean_k is not None else 'None'},\n"
            f"  s_fk_squared={'computed' if self.s_fk_squared is not None else 'None'}\n",
            f"  N_mean_k={'computed' if self.N_mean_k is not None else 'None'},\n"
            f"  s_Nk_squared={'computed' if self.s_Nk_squared is not None else 'None'}"
            f")"
        )

class GroundTruthDirichlet:
    """
    Object for the re-estimated Dirichlet distribution statistics.

    Parameters
    ----------
    f_mean_k : pandas.Series
        Estimated mean frequency for each mutation channel k.

    s_fk_squared : pandas.Series
        Estimated sample variance for each mutation channel k.

    Notes
    -----
    This class is intentionally lightweight. Its purpose is to group
    related results into a single object with explicit attribute names,
    rather than passing multiple pandas Series independently between
    functions.

    The two Series are expected to refer to the same mutation channels
    and to have compatible indices.
    """

    def __init__(self, f_mean_k: pd.Series, s_fk_squared: pd.Series):
        self.f_mean_k = f_mean_k
        self.s_fk_squared = s_fk_squared

class GroundTruthDirichletMle:
    """Object for storing sample frequencies and automatically deriving channel-level

    statistics for Dirichlet Maximum Likelihood Estimation (MLE).

    Parameters
    ----------
    f_kg : pandas.DataFrame
        Matrix of sample mutation frequencies. Rows correspond to mutation
        channels k (96 channels) and columns correspond to sample indices g (G
        samples).
    """

    def __init__(self, f_kg: pd.DataFrame):
        self.f_kg = f_kg
        self.f_mean_k = self.f_kg.mean(axis=1)
        self.s_fk_squared = self.f_kg.var(axis=1, ddof=1)
        self.log_f_kg = np.log(self.f_kg)
        self.mean_log_f_k = self.log_f_kg.mean(axis=1)

    def __repr__(self):
        return (
            f"GroundTruthDirichletMle("
            f"channels={self.f_kg.shape[0]}, "
            f"samples={self.f_kg.shape[1]})"
        )

class GroundTruthNegativeBinomialHierarchical:
    """
    Object for the re-estimated statistics from the hierarchical
    Negative Binomial model.

    Parameters
    ----------
    N_mean : float
        Estimated mean mutation burden across the samples.

    s_N_squared : float
        Estimated sample variance of the mutation burden.

    Notes
    -----
    This class represents the hierarchical Negative Binomial model,
    where mutation counts are modeled through a shared distribution
    of sample-level mutation burdens.

    The class provides a consistent object-based interface for storing
    the statistics required by the evaluation framework.
    """

    def __init__(self, N_mean: float, s_N_squared: float):
        self.N_mean = N_mean
        self.s_N_squared = s_N_squared


class GroundTruthNegativeBinomialNonHierarchical:
    """
    Object for the re-estimated statistics from the non-hierarchical
    Negative Binomial model.

    Parameters
    ----------
    N_mean_k : pandas.Series
        Estimated mean mutation count for each mutation channel k.

    s_Nk_squared : pandas.Series
        Estimated sample variance of the mutation count for each
        mutation channel k.

    Notes
    -----
    Unlike the hierarchical model, this representation stores
    channel-specific mean and variance estimates. The indices of
    N_mean_k and s_Nk_squared should correspond to the same
    mutation channels.
    """

    def __init__(self, N_mean_k: pd.Series, s_Nk_squared: pd.Series):
        self.N_mean_k = N_mean_k
        self.s_Nk_squared = s_Nk_squared


class GroundTruthMultinomial:
    """
    Object for the re-estimated statistics from the Multinomial model.

    Parameters
    ----------
    N_g : float
        Estimated total mutation count for group g.

    f_kg : pandas.Series
        Estimated mutation-frequency vector for group g, where
        each element corresponds to a mutation channel k.

    Notes
    -----
    The total count N_g and frequency vector f_kg jointly
    describe the Multinomial distribution for the group.
    """

    def __init__(self, N_g: float, f_kg: pd.Series):
        self.N_g = N_g
        self.f_kg = f_kg