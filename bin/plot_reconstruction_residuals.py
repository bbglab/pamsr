#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_reconstruction_residuals(
    csv_file: str,
    output_plot: str,
    color: str = "mediumpurple",
) -> None:
    """
    Computes and plots reconstruction residuals (N_with - N_without) across injected mutation levels
    to evaluate potential count inflation or model overfitting.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing `n_injected`, `n_reconstructed_with`, and `n_reconstructed_without`.
    output_plot : str
        File path where the output PNG plot will be saved.
    color : str, default="mediumpurple"
        Color for the mean residual line plot and standard deviation shading.

    Returns
    -------
    None
        Saves the output figure directly to `output_plot`.

    Raises
    ------
    FileNotFoundError
        If `csv_file` does not exist at the specified path.
    KeyError
        If mandatory columns `n_injected`, `n_reconstructed_with`, or `n_reconstructed_without` are missing.

    Notes
    -----
    - Calculates sample-level residual: $\\Delta N_{reconstructed} = N_{with} - N_{without}$.
    - Groups by `n_injected` to plot mean residual and $\\pm 1$ SD shaded bounds.
    - Includes a reference baseline at $y = 0$.

    Examples
    --------
    >>> plot_reconstruction_residuals(
    ...     csv_file="output/msigact/summary_stats.csv",
    ...     output_plot="output/msigact/reconstruction_residuals_5000.png",
    ...     color="mediumpurple",
    ... )
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Input CSV file not found: {csv_file}")

    df = pd.read_csv(csv_file)

    required_cols = ["n_injected", "n_reconstructed_with", "n_reconstructed_without"]
    for col in required_cols:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' missing from {csv_file}")

    df["n_injected"] = pd.to_numeric(df["n_injected"], errors="coerce")
    df = df.dropna(subset=["n_injected"])

    # Compute reconstruction residual per sample
    df["residual"] = df["n_reconstructed_with"] - df["n_reconstructed_without"]

    # Calculate summary statistics per n_injected group
    summary = (
        df.groupby("n_injected")["residual"]
        .agg(mean="mean", std="std")
        .reset_index()
    )
    summary["std"] = summary["std"].fillna(0)

    plt.figure(figsize=(10, 5))
    sns.set_theme(style="whitegrid")

    # Reference line at zero residual
    plt.axhline(
        0, color="gray", linestyle="--", linewidth=1.5, label="Zero Bias Line"
    )

    # Plot mean line and SD shading
    plt.plot(
        summary["n_injected"],
        summary["mean"],
        marker="o",
        color=color,
        linewidth=2,
        label=r"Mean Residual ($\Delta N$)",
    )
    plt.fill_between(
        summary["n_injected"],
        summary["mean"] - summary["std"],
        summary["mean"] + summary["std"],
        color=color,
        alpha=0.2,
        label=r"$\pm 1$ SD",
    )

    plt.title(
        r"Reconstruction Residuals ($\Delta N_{reconstructed} = N_{with} - N_{without}$)",
        fontsize=13,
    )
    plt.xlabel(r"Number of Injected Mutations ($n\_injected$)", fontsize=11)
    plt.ylabel(r"$\Delta N_{reconstructed}$ (Mutations)", fontsize=11)
    plt.legend(loc="best")
    plt.tight_layout()

    # Save output plot
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot reconstruction residuals across injected mutation levels."
    )
    parser.add_argument(
        "--csv_file",
        type=str,
        required=True,
        help="Path to summary statistics CSV file.",
    )
    parser.add_argument(
        "--output_plot",
        type=str,
        required=True,
        help="Path to save output plot PNG.",
    )
    parser.add_argument(
        "--color",
        type=str,
        default="mediumpurple",
        help="Color for mean line and SD shading (default: mediumpurple).",
    )

    args = parser.parse_args()

    plot_reconstruction_residuals(
        csv_file=args.csv_file,
        output_plot=args.output_plot,
        color=args.color,
    )