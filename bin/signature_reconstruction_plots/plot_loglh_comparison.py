#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_loglh_comparison(
    csv_file: str,
    output_plot: str,
    color_with: str = "dodgerblue",
    color_without: str = "darkorange",
    stat: str = "median",
) -> None:
    """
    Plots and compares log-likelihood values (`loglh_with` vs `loglh_without`)
    across varying numbers of injected mutations.

    Parameters
    ----------
    csv_file : str
        Path to the summary statistics CSV file.
    output_plot : str
        File path where the output PNG plot will be saved.
    color_with : str, default="dodgerblue"
        Color for the 'With Target' model line and shading.
    color_without : str, default="darkorange"
        Color for the 'Without Target' model line and shading.
    stat : {"median", "mean"}, default="median"
        Summary statistic computed per injection level ('median' with IQR shading or 'mean' with ±1 SD shading).

    Returns
    -------
    None
        Saves the output figure directly to `output_plot`.

    Raises
    ------
    FileNotFoundError
        If `csv_file` does not exist at the specified path.
    KeyError
        If mandatory columns `n_injected`, `loglh_with`, or `loglh_without` are missing.
    ValueError
        If `stat` is not one of 'median' or 'mean'.

    Notes
    -----
    - Aggregates metrics by `n_injected` groups and computes confidence bounds via percentile or standard deviation.

    Examples
    --------
    >>> plot_loglh_comparison(
    ...     csv_file="output/msigact/summary_stats.csv",
    ...     output_plot="output/msigact/comparison_log_likelihood_5000.png",
    ...     stat="median",
    ... )
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Input CSV file not found: {csv_file}")

    if stat not in ["median", "mean"]:
        raise ValueError("Parameter 'stat' must be either 'median' or 'mean'.")

    df = pd.read_csv(csv_file)

    for col in ["n_injected", "loglh_with", "loglh_without"]:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' missing from {csv_file}")

    df["n_injected"] = pd.to_numeric(df["n_injected"], errors="coerce")
    df = df.dropna(subset=["n_injected"])

    plt.figure(figsize=(10, 5))
    sns.set_theme(style="whitegrid")

    if stat == "median":
        agg_funcs = {
            "center": "median",
            "lower": lambda x: x.quantile(0.25),
            "upper": lambda x: x.quantile(0.75),
        }
        stat_label = "Median"
        shade_label = "IQR (25th-75th)"
    else:
        agg_funcs = {
            "center": "mean",
            "lower": lambda x: x.mean() - x.std(),
            "upper": lambda x: x.mean() + x.std(),
        }
        stat_label = "Mean"
        shade_label = "±1 SD"

    # Compute statistics for both models
    summary_with = (
        df.groupby("n_injected")["loglh_with"].agg(**agg_funcs).reset_index()
    )
    summary_without = (
        df.groupby("n_injected")["loglh_without"]
        .agg(**agg_funcs)
        .reset_index()
    )

    # Plot 'With Target' model
    plt.plot(
        summary_with["n_injected"],
        summary_with["center"],
        marker="o",
        linewidth=2,
        color=color_with,
        label=f"With Target ({stat_label})",
    )
    plt.fill_between(
        summary_with["n_injected"],
        summary_with["lower"].fillna(summary_with["center"]),
        summary_with["upper"].fillna(summary_with["center"]),
        color=color_with,
        alpha=0.2,
        label=f"With Target ({shade_label})",
    )

    # Plot 'Without Target' model
    plt.plot(
        summary_without["n_injected"],
        summary_without["center"],
        marker="s",
        linewidth=2,
        linestyle="--",
        color=color_without,
        label=f"Without Target ({stat_label})",
    )
    plt.fill_between(
        summary_without["n_injected"],
        summary_without["lower"].fillna(summary_without["center"]),
        summary_without["upper"].fillna(summary_without["center"]),
        color=color_without,
        alpha=0.2,
        label=f"Without Target ({shade_label})",
    )

    plt.title(
        f"Log-Likelihood Comparison: With vs. Without Target Signature ({stat_label})",
        fontsize=13,
    )
    plt.xlabel("Number of Injected Mutations ($n\\_injected$)", fontsize=11)
    plt.ylabel("Log-Likelihood ($\log L$)", fontsize=11)
    plt.legend(loc="best")
    plt.tight_layout()

    # Save figure
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot log-likelihood comparison between 'With' and 'Without' target models."
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
        "--color_with",
        type=str,
        default="dodgerblue",
        help="Color for 'With Target' model curve.",
    )
    parser.add_argument(
        "--color_without",
        type=str,
        default="darkorange",
        help="Color for 'Without Target' model curve.",
    )
    parser.add_argument(
        "--stat",
        type=str,
        choices=["median", "mean"],
        default="median",
        help="Aggregation statistic: 'median' or 'mean' (default: median).",
    )

    args = parser.parse_args()

    plot_loglh_comparison(
        csv_file=args.csv_file,
        output_plot=args.output_plot,
        color_with=args.color_with,
        color_without=args.color_without,
        stat=args.stat,
    )