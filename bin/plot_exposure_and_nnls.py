#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import nnls


def plot_exposure_and_nnls(
    cosmic_file: str,
    exposures_file: str,
    target_sig: str,
    bg_sigs: list,
    output_plot: str,
) -> None:
    """
    Computes Non-Negative Least Squares (NNLS) decomposition for a target SBS signature 
    against specified background signatures and plots exposure trends across injected mutations.

    Parameters
    ----------
    cosmic_file : str
        Path to the CSV file containing reference COSMIC signature profiles.
    exposures_file : str
        Path to the CSV file containing target signature exposures and mutation counts.
    target_sig : str
        Target SBS signature identifier (e.g., "SBS31").
    bg_sigs : list of str
        List of background signature identifiers (e.g., ["SBS1", "SBS5", "SBS40"]).
    output_plot : str
        File path where the output PNG plot will be saved.

    Returns
    -------
    None
        Saves the output figure directly to `output_plot` and prints weights and residual to stdout.

    Raises
    ------
    FileNotFoundError
        If `cosmic_file` or `exposures_file` does not exist at the specified path.
    KeyError
        If `target_sig` or any column in `bg_sigs` is missing from the input dataframes.

    Notes
    -----
    - The first column of the COSMIC file is used as the row index for mutation contexts.
    - Profiles are converted to numeric matrices and solved via scipy `nnls`.
    - Includes a $y = x$ reference line in the plot.
    - Displays the computed residual strictly within the plot legend without rendering an additional trendline.

    Examples
    --------
    >>> plot_exposure_and_nnls(
    ...     cosmic_file="input/cosmic_v3_6_sbs.csv",
    ...     exposures_file="output/exposures_target.csv",
    ...     target_sig="SBS31",
    ...     bg_sigs=["SBS1", "SBS5", "SBS40"],
    ...     output_plot="output/exposure_plot.png",
    ... )
    """
    if not os.path.exists(cosmic_file):
        raise FileNotFoundError(f"COSMIC file not found: {cosmic_file}")
    if not os.path.exists(exposures_file):
        raise FileNotFoundError(f"Exposures file not found: {exposures_file}")

    # Load data
    df_cosmic = pd.read_csv(cosmic_file, sep="\t")
    context_col = df_cosmic.columns[0]
    df_cosmic.set_index(context_col, inplace=True)

    df_exp = pd.read_csv(exposures_file)

    # Validate signature columns
    missing_bg = [sig for sig in bg_sigs if sig not in df_cosmic.columns]
    if missing_bg:
        raise KeyError(
            f"Background signature(s) {missing_bg} not found in COSMIC file."
        )

    if target_sig not in df_cosmic.columns:
        raise KeyError(
            f"Target signature '{target_sig}' not found in COSMIC file."
        )

    if target_sig not in df_exp.columns:
        raise KeyError(
            f"Target column '{target_sig}' not found in exposures file."
        )

    # NNLS Calculation
    b_matrix = (
        df_cosmic[bg_sigs]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .to_numpy()
    )
    t_vector = (
        df_cosmic[target_sig]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .to_numpy()
    )

    weights, residual = nnls(A=b_matrix, b=t_vector)

    print(f"NNLS Weights: {weights}")
    print(f"NNLS Residual: {residual}")

    # Plot Setup
    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid")

    ax = sns.lineplot(
        data=df_exp,
        x="n_injected",
        y=target_sig,
        errorbar="sd",
        color="#006AFF",
        marker="o",
        label=f"{target_sig} Exposure",
    )

    # Baseline y = x line
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    lower = min(xmin, ymin)
    upper = max(xmax, ymax)

    ax.plot(
        [lower, upper],
        [lower, upper],
        linestyle="--",
        color="black",
        linewidth=1.5,
        label="y = x",
    )

    # Include NNLS residual strictly in the legend without plotting a line
    ax.plot([], [], " ", label=f"NNLS Residual: {residual:.4f}")
    ax.legend(title="Legend", loc="best")

    plt.title(
        f"Exposure of {target_sig} by Number of Injected Mutations",
        fontsize=14,
    )
    plt.xlabel("n_injected", fontsize=12)
    plt.ylabel(f"{target_sig} Value (Mean ± SD)", fontsize=12)

    plt.tight_layout()

    # Save output plot
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot signature exposures and NNLS residual using named signature columns."
    )
    parser.add_argument(
        "--cosmic_file",
        type=str,
        required=True,
        help="Path to the COSMIC signatures CSV file.",
    )
    parser.add_argument(
        "--exposures_file",
        type=str,
        required=True,
        help="Path to the exposures CSV file.",
    )
    parser.add_argument(
        "--target_sig",
        type=str,
        required=True,
        help="Target signature name (e.g., SBS31).",
    )
    parser.add_argument(
        "--bg_sigs",
        type=str,
        nargs="+",
        required=True,
        help="List of background signatures (e.g., SBS1 SBS2 SBS6).",
    )
    parser.add_argument(
        "--output_plot",
        type=str,
        required=True,
        help="Path to save output plot PNG.",
    )

    args = parser.parse_args()

    plot_exposure_and_nnls(
        cosmic_file=args.cosmic_file,
        exposures_file=args.exposures_file,
        target_sig=args.target_sig,
        bg_sigs=args.bg_sigs,
        output_plot=args.output_plot,
    )