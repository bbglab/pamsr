import numpy as np
import pandas as pd

TRIMUT_TO_TRI={'A[C>A]A': 'ACA',
 'A[C>G]A': 'ACA',
 'A[C>T]A': 'ACA',
 'A[C>A]C': 'ACC',
 'A[C>G]C': 'ACC',
 'A[C>T]C': 'ACC',
 'A[C>A]G': 'ACG',
 'A[C>G]G': 'ACG',
 'A[C>T]G': 'ACG',
 'A[C>A]T': 'ACT',
 'A[C>G]T': 'ACT',
 'A[C>T]T': 'ACT',
 'A[T>A]A': 'ATA',
 'A[T>C]A': 'ATA',
 'A[T>G]A': 'ATA',
 'A[T>A]C': 'ATC',
 'A[T>C]C': 'ATC',
 'A[T>G]C': 'ATC',
 'A[T>A]G': 'ATG',
 'A[T>C]G': 'ATG',
 'A[T>G]G': 'ATG',
 'A[T>A]T': 'ATT',
 'A[T>C]T': 'ATT',
 'A[T>G]T': 'ATT',
 'C[C>A]A': 'CCA',
 'C[C>G]A': 'CCA',
 'C[C>T]A': 'CCA',
 'C[C>A]C': 'CCC',
 'C[C>G]C': 'CCC',
 'C[C>T]C': 'CCC',
 'C[C>A]G': 'CCG',
 'C[C>G]G': 'CCG',
 'C[C>T]G': 'CCG',
 'C[C>A]T': 'CCT',
 'C[C>G]T': 'CCT',
 'C[C>T]T': 'CCT',
 'C[T>A]A': 'CTA',
 'C[T>C]A': 'CTA',
 'C[T>G]A': 'CTA',
 'C[T>A]C': 'CTC',
 'C[T>C]C': 'CTC',
 'C[T>G]C': 'CTC',
 'C[T>A]G': 'CTG',
 'C[T>C]G': 'CTG',
 'C[T>G]G': 'CTG',
 'C[T>A]T': 'CTT',
 'C[T>C]T': 'CTT',
 'C[T>G]T': 'CTT',
 'G[C>A]A': 'GCA',
 'G[C>G]A': 'GCA',
 'G[C>T]A': 'GCA',
 'G[C>A]C': 'GCC',
 'G[C>G]C': 'GCC',
 'G[C>T]C': 'GCC',
 'G[C>A]G': 'GCG',
 'G[C>G]G': 'GCG',
 'G[C>T]G': 'GCG',
 'G[C>A]T': 'GCT',
 'G[C>G]T': 'GCT',
 'G[C>T]T': 'GCT',
 'G[T>A]A': 'GTA',
 'G[T>C]A': 'GTA',
 'G[T>G]A': 'GTA',
 'G[T>A]C': 'GTC',
 'G[T>C]C': 'GTC',
 'G[T>G]C': 'GTC',
 'G[T>A]G': 'GTG',
 'G[T>C]G': 'GTG',
 'G[T>G]G': 'GTG',
 'G[T>A]T': 'GTT',
 'G[T>C]T': 'GTT',
 'G[T>G]T': 'GTT',
 'T[C>A]A': 'TCA',
 'T[C>G]A': 'TCA',
 'T[C>T]A': 'TCA',
 'T[C>A]C': 'TCC',
 'T[C>G]C': 'TCC',
 'T[C>T]C': 'TCC',
 'T[C>A]G': 'TCG',
 'T[C>G]G': 'TCG',
 'T[C>T]G': 'TCG',
 'T[C>A]T': 'TCT',
 'T[C>G]T': 'TCT',
 'T[C>T]T': 'TCT',
 'T[T>A]A': 'TTA',
 'T[T>C]A': 'TTA',
 'T[T>G]A': 'TTA',
 'T[T>A]C': 'TTC',
 'T[T>C]C': 'TTC',
 'T[T>G]C': 'TTC',
 'T[T>A]G': 'TTG',
 'T[T>C]G': 'TTG',
 'T[T>G]G': 'TTG',
 'T[T>A]T': 'TTT',
 'T[T>C]T': 'TTT',
 'T[T>G]T': 'TTT'}

GCHR37_CONTEXT={'ACA': 0.04033784990481213,
 'ACC': 0.023253615661324983,
 'ACG': 0.005021484955974482,
 'ACT': 0.03217912683331554,
 'ATA': 0.0412399330162748,
 'ATC': 0.026699871287136336,
 'ATG': 0.036732020382732364,
 'ATT': 0.0498693509297776,
 'CCA': 0.03688223848637286,
 'CCC': 0.026284566095707965,
 'CCG': 0.00551552607744207,
 'CCT': 0.03551064955287565,
 'CTA': 0.02579312514648204,
 'CTC': 0.03365706991273405,
 'CTG': 0.04051762874504749,
 'CTT': 0.03991192818878235,
 'GCA': 0.028800666334767897,
 'GCC': 0.023789566326101134,
 'GCG': 0.004755495346407583,
 'GCT': 0.0279590675743442,
 'GTA': 0.022690906625059314,
 'GTC': 0.018887943258674038,
 'GTG': 0.03005352468027329,
 'GTT': 0.02915970279142049,
 'TCA': 0.0391755019004158,
 'TCC': 0.030865231429635178,
 'TCG': 0.004410821229553994,
 'TCT': 0.04423091572843269,
 'TTA': 0.04161022370962304,
 'TTC': 0.03943758547967859,
 'TTG': 0.037893085267018047,
 'TTT': 0.076873777141802}

GCHR38_CONTEXT={'ACA': 0.04050331474286235,
 'ACC': 0.02299638005142609,
 'ACG': 0.005177692139446636,
 'ACT': 0.03230325620576059,
 'ATA': 0.0410345979390114,
 'ATC': 0.026766974486020972,
 'ATG': 0.0365675394955667,
 'ATT': 0.049672832908862284,
 'CCA': 0.036630583120720164,
 'CCC': 0.025923477399463137,
 'CCG': 0.005526945658854522,
 'CCT': 0.03525962926004756,
 'CTA': 0.025747018693359606,
 'CTC': 0.03370328474071515,
 'CTG': 0.04044547103103616,
 'CTT': 0.04009979894724738,
 'GCA': 0.02881111694243897,
 'GCC': 0.023471044467824956,
 'GCG': 0.004794732373654194,
 'GCT': 0.027826534257835513,
 'GTA': 0.022520149630070003,
 'GTC': 0.018787228927003696,
 'GTG': 0.03014663057082029,
 'GTT': 0.029526626856226548,
 'TCA': 0.03934099901391775,
 'TCC': 0.030949722276210286,
 'TCG': 0.004498105275405602,
 'TCT': 0.0446061070084102,
 'TTA': 0.04103538332661506,
 'TTC': 0.04013744848679336,
 'TTG': 0.03812637647056972,
 'TTT': 0.07706299729580315}

EQUIPROBABLE_CONTEXT={'ACA': 0.010416666666666666,
 'ACC': 0.010416666666666666,
 'ACG': 0.010416666666666666,
 'ACT': 0.010416666666666666,
 'ATA': 0.010416666666666666,
 'ATC': 0.010416666666666666,
 'ATG': 0.010416666666666666,
 'ATT': 0.010416666666666666,
 'CCA': 0.010416666666666666,
 'CCC': 0.010416666666666666,
 'CCG': 0.010416666666666666,
 'CCT': 0.010416666666666666,
 'CTA': 0.010416666666666666,
 'CTC': 0.010416666666666666,
 'CTG': 0.010416666666666666,
 'CTT': 0.010416666666666666,
 'GCA': 0.010416666666666666,
 'GCC': 0.010416666666666666,
 'GCG': 0.010416666666666666,
 'GCT': 0.010416666666666666,
 'GTA': 0.010416666666666666,
 'GTC': 0.010416666666666666,
 'GTG': 0.010416666666666666,
 'GTT': 0.010416666666666666,
 'TCA': 0.010416666666666666,
 'TCC': 0.010416666666666666,
 'TCG': 0.010416666666666666,
 'TCT': 0.010416666666666666,
 'TTA': 0.010416666666666666,
 'TTC': 0.010416666666666666,
 'TTG': 0.010416666666666666,
 'TTT': 0.010416666666666666}

def modify_trinucleotide_context(
    df: pd.DataFrame,
    mutation_types: list,
    orig_context: dict,
    target_context: dict,
) -> pd.DataFrame:
    """Normalizes a 96-mutation count matrix from an original trinucleotide context

    C to a target context C_tilde using dictionary mapping.
    """
    res_df = df.copy()

    # Direct lookup from 96 mutation channels to 32 trinucleotide contexts
    trinucleotides = pd.Series(mutation_types).map(TRIMUT_TO_TRI)

    # Extract original and target background context proportions
    c_orig = trinucleotides.map(orig_context).astype(float).values
    c_target = trinucleotides.map(target_context).astype(float).values

    # Perform matrix context normalization directly on sample columns
    for col in res_df.columns:
        s = res_df[col].values.astype(float)
        # Remove original background context bias
            # out argument indicates to write zero in all positions
            # where argument indicates where to apply the actual calculation
        s_prime = np.divide(s, c_orig, out=np.zeros_like(s), where=c_orig != 0)
        # Scale by target background proportions
        s_target_unnorm = s_prime * c_target
        # Re-normalize so sample proportions sum to 1
        total_weight = np.sum(s_target_unnorm)
        if total_weight > 0:
            res_df[col] = s_target_unnorm / total_weight
        else:
            res_df[col] = 0.0
            
    return res_df

import pandas as pd
import numpy as np

def reconstruct_matrix_without_signature(
    df_prob: pd.DataFrame,
    df_act: pd.DataFrame,
    target_signature: str,
    mutation_col: str = "MutationType",
    sample_col: str = "Samples"
) -> pd.DataFrame:
    """
    Brief description
    -----------------
    Reconstructs a mutational matrix (M = P x E^T) after zeroing out the 
    exposure/activity of a specified target signature.

    Parameters
    ----------
    df_prob : pd.DataFrame
        DataFrame of signature probabilities (P) of shape (MutationTypes x Signatures).
    df_act : pd.DataFrame
        DataFrame of signature activities/exposures (E) of shape (Samples x Signatures).
    target_signature : str
        The name of the signature to eliminate by setting its exposures to zero.
    mutation_col : str, default="MutationType"
        Column name for mutation types in `df_prob` if not already the index.
    sample_col : str, default="Samples"
        Column name for sample identifiers in `df_act` if not already the index.

    Returns
    -------
    pd.DataFrame
        Reconstructed mutational matrix (M) of shape (MutationTypes x Samples).

    Raises
    ------
    ValueError
        If no overlapping signatures exist between `df_prob` and `df_act`.
    ValueError
        If `target_signature` is not found within the common signature columns.

    Notes
    -----
    - Leading and trailing whitespaces are automatically stripped from signature column headers.
    - Signatures present in only one of the input matrices are dropped during column alignment.
    - Zeroing out the exposure column in E ensures that the contribution of `target_signature` 
      is removed while maintaining original matrix dimensions and sample order.

    Examples
    --------
    >>> import pandas as pd
    >>> df_p = pd.DataFrame({
    ...     'MutationType': ['C>A', 'C>G'],
    ...     'SBS1': [0.8, 0.2],
    ...     'SBS5': [0.5, 0.5]
    ... })
    >>> df_e = pd.DataFrame({
    ...     'Samples': ['Sample_A', 'Sample_B'],
    ...     'SBS1': [100, 200],
    ...     'SBS5': [50, 75]
    ... })
    >>> reconstructed_M = reconstruct_matrix_without_signature(
    ...     df_prob=df_p,
    ...     df_act=df_e,
    ...     target_signature='SBS1'
    ... )
    >>> print(reconstructed_M)
    Samples       Sample_A  Sample_B
    MutationType                  
    C>A               25.0     37.5
    C>G               25.0     37.5
    """
    # 1. Prepare indices
    P = df_prob.set_index(mutation_col) if mutation_col in df_prob.columns else df_prob.copy()
    E = df_act.set_index(sample_col) if sample_col in df_act.columns else df_act.copy()

    # 2. Clean whitespace from columns
    P.columns = P.columns.str.strip()
    E.columns = E.columns.str.strip()
    target_sig_clean = target_signature.strip()

    # 3. Align on shared signatures
    common_sigs = P.columns.intersection(E.columns)
    if common_sigs.empty:
        raise ValueError("No common signature columns found between df_prob and df_act.")

    if target_sig_clean not in common_sigs:
        raise ValueError(
            f"Target signature '{target_sig_clean}' not found in shared signatures: {list(common_sigs)}"
        )

    P_aligned = P[common_sigs].copy()
    E_aligned = E[common_sigs].copy()

    # 4. Zero out the target signature exposure in E
    E_aligned[target_sig_clean] = 0.0

    # 5. Perform matrix multiplication (M = P x E^T)
    M_values = P_aligned.values @ E_aligned.T.values

    # 6. Construct output DataFrame
    M_df = pd.DataFrame(M_values, index=P_aligned.index, columns=E_aligned.index)
    M_df.index.name = P_aligned.index.name or "MutationType"

    # --- Convert floating mutation counts to integers ---
    return M_df.round().astype(int)

reconstruct_matrix_without_signature(df_prob,df_act,"SBS1")