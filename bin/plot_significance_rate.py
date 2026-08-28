#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_significance_rate(
    csv_file: str,
    output_plot: str,
    alpha: float = 0.05,
    color: str = "crimson",
) -> None:
    """
    Calculates and plots the percentage of samples meeting a chi-square p-value
    significance threshold against the number of injected mutations.

    Parameters
    ----------
    csv_file : str
        Path to the summary statistics CSV file containing `n_injected` and `chisq_p`.
    output_plot : str
        File path where the output PNG plot will be saved.
    alpha : float, default=0.05
        Significance threshold cutoff for p-values (`chisq_p < alpha`).
    color : str, default="crimson"
        Color for the line plot and data points.

    Returns
    -------
    None
        Saves the output figure directly to `output_plot`.

    Raises
    ------
    FileNotFoundError
        If `csv_file` does not exist at the specified path.
    KeyError
        If mandatory columns `n_injected` or `chisq_p` are missing from `csv_file`.

    Notes
    -----
    - Computes sample percentage meeting significance via `(chisq_p < alpha).mean() * 100`.
    - Automatically creates missing parent directories for the output path.

    Examples
    --------
    >>> plot_significance_rate(
    ...     csv_file="output/msigact/summary_stats.csv",
    ...     output_plot="output/msigact/percentage_p_value_5000.png",
    ...     alpha=0.05,
    ...     color="crimson",
    ... )
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Input CSV file not found: {csv_file}")

    df = pd.read_csv(csv_file)

    for col in ["n_injected", "chisq_p"]:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' missing from {csv_file}")

    df["n_injected"] = pd.to_numeric(df["n_injected"], errors="coerce")
    df = df.dropna(subset=["n_injected"])

    # Calculate percentage of samples meeting significance threshold
    summary_df = (
        df.groupby("n_injected")["chisq_p"]
        .apply(lambda p_vals: (p_vals < alpha).mean() * 100)
        .reset_index(name="pct_significant")
    )

    plt.figure(figsize=(9, 5))
    sns.set_theme(style="whitegrid")

    plt.plot(
        summary_df["n_injected"],
        summary_df["pct_significant"],
        marker="o",
        linewidth=2,
        color=color,
        label=f"p < {alpha}",
    )

    plt.title(
        f"Percentage of Significant Samples ($p < {alpha}$) vs Injected Mutations",
        fontsize=13,
    )
    plt.xlabel("Number of Injected Mutations ($n\\_injected$)", fontsize=11)
    plt.ylabel(f"% of Samples ($p < {alpha}$)", fontsize=11)
    plt.ylim(-5, 105)
    plt.tight_layout()

    # Save figure
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot percentage of significant p-value samples across injected mutations."
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
        "--alpha",
        type=float,
        default=0.05,
        help="Significance threshold for p-value (default: 0.05).",
    )
    parser.add_argument(
        "--color",
        type=str,
        default="crimson",
        help="Plot line and marker color (default: crimson).",
    )

    args = parser.parse_args()

    plot_significance_rate(
        csv_file=args.csv_file,
        output_plot=args.output_plot,
        alpha=args.alpha,
        color=args.color,
    )