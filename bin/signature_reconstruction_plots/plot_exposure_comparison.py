#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_exposure_comparison(
    file_with: str,
    file_without: str,
    target_sig: str,
    output_plot: str,
    color_with: str = "blue",
    color_without: str = "red",
) -> None:
    """
    Plots a line plot comparing the exposure of a target signature across injected mutations
    between datasets generated with and without the target signature.

    Parameters
    ----------
    file_with : str
        Path to the CSV file containing exposures with the target signature present.
    file_without : str
        Path to the CSV file containing exposures without the target signature present.
    target_sig : str
        Target SBS signature identifier (e.g., "SBS3") to compare.
    output_plot : str
        File path where the output PNG plot will be saved.
    color_with : str, default="blue"
        Color for the "With Target" dataset line plot.
    color_without : str, default="red"
        Color for the "Without Target" dataset line plot.

    Returns
    -------
    None
        Saves the output figure directly to `output_plot`.

    Raises
    ------
    FileNotFoundError
        If `file_with` or `file_without` does not exist at the specified path.
    KeyError
        If `target_sig` or 'n_injected' columns are missing from either input file.

    Notes
    -----
    - Uses Seaborn lineplot with standard deviation error bars (`errorbar='sd'`).
    - Creates destination directory automatically if it does not already exist.

    Examples
    --------
    >>> plot_exposure_comparison(
    ...     file_with="output/msigact/exposures_with_target_5000.csv",
    ...     file_without="output/msigact/exposures_without_target_5000.csv",
    ...     target_sig="SBS3",
    ...     output_plot="output/msigact/comparing_sign_5000.png",
    ... )
    """
    if not os.path.exists(file_with):
        raise FileNotFoundError(f"Input file not found: {file_with}")
    if not os.path.exists(file_without):
        raise FileNotFoundError(f"Input file not found: {file_without}")

    # Load data
    df_with = pd.read_csv(file_with)
    df_without = pd.read_csv(file_without)

    # Validate required columns
    for df_name, df in [("With Target", df_with), ("Without Target", df_without)]:
        if "n_injected" not in df.columns:
            raise KeyError(
                f"Column 'n_injected' missing in {df_name} file."
            )
        if target_sig not in df.columns:
            raise KeyError(
                f"Target column '{target_sig}' missing in {df_name} file."
            )

    # Plot setup
    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid")

    # Plot lineplots
    sns.lineplot(
        data=df_with,
        x="n_injected",
        y=target_sig,
        errorbar="sd",
        color=color_with,
        marker="o",
        label="With Target",
    )

    sns.lineplot(
        data=df_without,
        x="n_injected",
        y=target_sig,
        errorbar="sd",
        color=color_without,
        marker="o",
        label="Without Target",
    )

    # Format plot
    plt.title(
        f"Comparison of {target_sig} Exposure (Mean ± SD)", fontsize=14
    )
    plt.xlabel("Number of Injected Mutations (n_injected)", fontsize=12)
    plt.ylabel(f"{target_sig} Value", fontsize=12)
    plt.legend(title="Dataset", loc="upper left")

    plt.tight_layout()

    # Save output plot
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare target signature exposures with and without target dataset."
    )
    parser.add_argument(
        "--file_with",
        type=str,
        required=True,
        help="Path to exposures CSV file with target signature.",
    )
    parser.add_argument(
        "--file_without",
        type=str,
        required=True,
        help="Path to exposures CSV file without target signature.",
    )
    parser.add_argument(
        "--target_sig",
        type=str,
        required=True,
        help="Target signature name to compare (e.g., SBS3).",
    )
    parser.add_argument(
        "--output_plot",
        type=str,
        required=True,
        help="Path to save output comparison plot PNG.",
    )
    parser.add_argument(
        "--color_with",
        type=str,
        default="blue",
        help="Color for the 'With Target' dataset curve.",
    )
    parser.add_argument(
        "--color_without",
        type=str,
        default="red",
        help="Color for the 'Without Target' dataset curve.",
    )

    args = parser.parse_args()

    plot_exposure_comparison(
        file_with=args.file_with,
        file_without=args.file_without,
        target_sig=args.target_sig,
        output_plot=args.output_plot,
        color_with=args.color_with,
        color_without=args.color_without,
    )


#     process COMPARE_EXPOSURES {
#     input:
#     path file_with_csv
#     path file_without_csv
#     val  target_signature

#     output:
#     path "*.png"

#     script:
#     """
#     compare_exposures.py \\
#         --file_with ${file_with_csv} \\
#         --file_without ${file_without_csv} \\
#         --target_sig ${target_signature} \\
#         --output_plot "comparison_${target_signature}.png"
#     """
# }