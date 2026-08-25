import numpy as np

from scipy.stats import dirichlet
from scipy.special import digamma, polygamma
from scipy.optimize import root_scalar

def estimate_w_mom(group):
    """
    Estimate the Dirichlet concentration weights (w_k) using the
    method of moments.

    Parameters
    ----------
    group : MutationGroup
        MutationGroup object with estimators already computed.

    Returns
    -------
    numpy.ndarray
        Estimated w_k for each mutation channel.
    """
    # Ensure estimators (f_kg, f_mean_k, s_fk_squared) are computed
    if group.f_mean_k is None:
        AttributeError("The f_mean_k has not be yet computed.")
        
    if group.s_fk_squared is None:
        AttributeError("The s_fk_squared has not be yet computed.")

    # Mean and Variance values are extracted
    mean = group.f_mean_k.values
    var = group.s_fk_squared.values

    valid = (var > 0) & (mean > 0) & (mean < 1)

    if not np.any(valid):
        raise ValueError("No valid channels with non-zero variance found. Try using the MLE mode for the computing the concentration parameter.")

    w = (mean * (1 - mean) / var) - 1

    return w

def estimate_w_mom_lse(group):
    """Estimates the optimal concentration weight W by minimizing the squared

    difference between theoretical MoM channel variance phi_i(W) and empirical
    channel variance S_i.

    Parameters:
    -----------
    group : MutationGroup
        Instance containing sample frequency matrix and computed estimators.

    Returns:
    --------
    W_opt : float
        Optimal least-squares concentration weight.
    """
    # Ensure estimators (f_kg, f_mean_k, s_fk_squared) are computed
    if group.f_mean_k is None:
        AttributeError("The f_mean_k has not be yet computed.")

    if group.s_fk_squared is None:
        AttributeError("The s_fk_squared has not be yet computed.")


    # Empirical expectations e_i and empirical variance S_i per channel
    e = group.f_mean_k.to_numpy()
    e = e / np.sum(e)  # Normalize to ensure sum(e) == 1
    S = group.s_fk_squared.to_numpy()

    # Mean-variance numerator factor: v_i = e_i * (1 - e_i)
    v = e * (1 - e)

    # Closed-form Least Squares solution:
    # Theoretical variance phi_i(W) = v_i / (W + 1) = c * v_i where c = 1 / (W + 1)
    # Minimizing sum((c * v_i - S_i)^2) w.r.t c yields c = sum(v_i * S_i) / sum(v_i^2)
    sum_v_squared = np.sum(v**2)
    sum_v_S = np.sum(v * S)

    if sum_v_S <= 0:
        return 1.0

    # Optimal scale factor c = 1 / (W + 1)
    c_opt = sum_v_S / sum_v_squared

    # Solve for W = (1 / c) - 1
    W_opt = (1.0 / c_opt) - 1.0

    return max(float(W_opt), 1.0)

def estimate_w_mle(group, eps=1e-12):
    """Finds the maximum likelihood scalar weight W given a MutationGroup instance.

    Parameters:
    -----------
    group : MutationGroup
        Instance containing sample frequency matrix and computed estimators.
    eps : float, default=1e-12
        Small constant to avoid log(0).

    Returns:
    --------
    W_opt : float
        Optimal MLE concentration weight.
    """
    if group.f_kg is None:
            raise AttributeError("The f_kg estimator has not been computed.")

    if group.f_mean_k is None:
        raise AttributeError("The f_mean_k estimator has not been computed.")

    if group.s_fk_squared is None:
        raise AttributeError("The s_fk_squared estimator has not been computed.")

    f_matrix = np.clip(group.f_kg.to_numpy(), eps, 1.0)

    e = group.f_mean_k.to_numpy()
    e = e / np.sum(e)

    # Ensure mean_log_f matches shape of e (K channels)
    if f_matrix.shape[0] == len(e):
        mean_log_f = np.mean(np.log(f_matrix), axis=1)
    else:
        mean_log_f = np.mean(np.log(f_matrix), axis=0)

    constant_term = np.sum(e * mean_log_f)

    # Objective and derivative in W-space
    def h(W):
        return digamma(W) - np.sum(e * digamma(W * e)) + constant_term

    def h_prime(W):
        return polygamma(1, W) - np.sum((e**2) * polygamma(1, W * e))

    # Objective and derivative in theta-space (theta = ln W)
    def g(theta):
        W = np.exp(theta)
        return h(W)

    def g_prime(theta):
        W = np.exp(theta)
        return W * h_prime(W)

    # Method of Moments initial guess
    var_f = group.s_fk_squared.to_numpy()
    mom_weights = (e * (1.0 - e)) / np.maximum(var_f, 1e-8) - 1.0
    W_init = np.maximum(np.median(mom_weights), 1.0)
    theta_init = np.log(W_init)

    sol = root_scalar(g, fprime=g_prime, x0=theta_init, method="newton")

    if not sol.converged:
        raise RuntimeError("MLE optimization for W failed to converge.")

    W_opt = float(np.exp(sol.root))
    return W_opt

def compute_error_variance_w_mom_lse(W, group):
    """Evaluates the Error Variance function EV(W) = sum_i (phi_i(W) - S_i)^2"""
    e = group.f_mean_k.to_numpy()
    e = e / np.sum(e)
    S = group.s_fk_squared.to_numpy()

    # Theoretical variance phi_i(W) under Method of Moments
    phi = (e * (1 - e)) / (W + 1)

    # Error variance loss function EV(W)
    return float(np.sum((phi - S) ** 2))

def estimate_theta(group,mode:str="mom"):
    """
    Estimate the Dirichlet theta parameters.

    The method-of-moments estimator is

        theta = W * mean_k

    Parameters
    ----------
    group : MutationGroup
    
    mode : str
        String defining if the weight parameter should be estimated by the method of moments or via maximum likelihood estimation.

        Possible values are:
            - mom
            - mle
            - mom_lse
    
    Returns
    -------
    numpy.ndarray
        Theta parameter for each mutation channel.
    """
    mean = group.f_mean_k.values

    if mode == "mom":
        w = estimate_w_mom(group)

    elif mode =="mom_lse":
        w = estimate_w_mom_lse(group)

    elif mode =="mle":
        w = estimate_w_mle(group)

    theta = w * mean
    
    return theta

def synthetic_frequencies(group, mode:str="mom", size:int=1, random_state:int=None):
    """
    Sample mutation frequency vectors from the estimated Dirichlet
    distribution.

    Parameters
    ----------
    group : MutationGroup

    size : int, default=1
        Number of synthetic frequency vectors to generate.

    random_state : int or None
        Random seed for reproducibility.

    Returns
    -------
    numpy.ndarray
        Samples from the Dirichlet distribution.
        Shape = (size, 96)
    """

    alpha = estimate_theta(group,mode)

    synth_freqs = dirichlet.rvs(
        alpha=alpha,
        size=size,
        random_state=random_state
    )

    return synth_freqs