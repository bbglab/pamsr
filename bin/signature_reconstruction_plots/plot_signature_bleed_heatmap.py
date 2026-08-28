#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_signature_bleed_heatmap(
    file_with: str,
    file_without: str,
    target_signature: str,
    output_plot: str,
    top_n: int = 15,
    cmap: str = "vlag",
) -> None:
    """
    Generates a heatmap illustrating signature crosstalk (delta exposure) across non-target
    signatures between datasets generated with and without a target signature.

    Parameters
    ----------
    file_with : str
        Path to the CSV file containing exposures with the target signature present.
    file_without : str
        Path to the CSV file containing exposures without the target signature present.
    target_signature : str
        Target SBS signature column identifier (e.g., "SBS3") to exclude from heatmap rows.
    output_plot : str
        File path where the output PNG heatmap will be saved.
    top_n : int, default=15
        Number of top impacted non-target signatures to display based on maximum absolute shift.
    cmap : str, default="vlag"
        Diverging color map used to represent negative and positive exposure shifts.

    Returns
    -------
    None
        Saves the output figure directly to `output_plot`.

    Raises
    ------
    FileNotFoundError
        If `file_with` or `file_without` does not exist at the specified path.
    KeyError
        If mandatory columns `sample` or `n_injected` are missing from input files.
    ValueError
        If no common non-target signatures are found between both input files.

    Notes
    -----
    - Subtraction (`df_with` minus `df_without`) is aligned strictly by sample identifier.
    - Symmetric color scaling (`vmin=-bound`, `vmax=bound`) centers the colorbar at zero shift.

    Examples
    --------
    >>> plot_signature_bleed_heatmap(
    ...     file_with="output/msigact/exposures_with_target_5000.csv",
    ...     file_without="output/msigact/exposures_without_target_5000.csv",
    ...     target_signature="SBS3",
    ...     output_plot="output/msigact/heatmap_delta_exposures_top15_5000.png",
    ...     top_n=15,
    ... )
    """
    if not os.path.exists(file_with):
        raise FileNotFoundError(f"Input file not found: {file_with}")
    if not os.path.exists(file_without):
        raise FileNotFoundError(f"Input file not found: {file_without}")

    # Load datasets
    df_with = pd.read_csv(file_with)
    df_without = pd.read_csv(file_without)

    meta_cols = ["sample", "n_injected"]
    for col in meta_cols:
        if col not in df_with.columns or col not in df_without.columns:
            raise KeyError(
                f"Mandatory metadata column '{col}' missing from input datasets."
            )

    # Identify common non-target signature columns
    common_non_target_sigs = [
        col
        for col in df_with.columns
        if col in df_without.columns
        and col not in meta_cols
        and col != target_signature
    ]

    if not common_non_target_sigs:
        raise ValueError(
            "No common non-target signature columns found between input files."
        )

    # Align datasets by sample for accurate subtraction
    df_with = df_with.sort_values("sample").reset_index(drop=True)
    df_without = df_without.sort_values("sample").reset_index(drop=True)

    # Compute Delta Exposure for common non-target signatures
    df_diff = (
        df_with[common_non_target_sigs] - df_without[common_non_target_sigs]
    )
    df_diff["n_injected"] = df_with["n_injected"]

    # Calculate mean exposure change grouped by n_injected
    mean_diff = df_diff.groupby("n_injected").mean()

    # Select top N signatures by maximum absolute shift
    top_sigs = mean_diff.abs().max().nlargest(top_n).index
    heatmap_data = mean_diff[top_sigs].T

    # Plot Heatmap centered at 0
    plt.figure(figsize=(11, 7))
    sns.set_theme(style="white")

    bound = max(abs(heatmap_data.values.min()), abs(heatmap_data.values.max()))

    sns.heatmap(
        heatmap_data,
        cmap=cmap,
        center=0,
        vmin=-bound,
        vmax=bound,
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": r"Mean $\Delta$ Exposure (With - Without)"},
    )

    plt.title(
        f"Signature Bleed / Crosstalk Heatmap (Target: {target_signature})\n"
        f"Top {top_n} Non-Target Signatures by Absolute Exposure Shift",
        fontsize=13,
    )
    plt.xlabel(r"Number of Injected Mutations ($n\_injected$)", fontsize=11)
    plt.ylabel("Non-Target SBS Signatures", fontsize=11)
    plt.tight_layout()

    # Save output plot
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate signature bleed/crosstalk heatmap for non-target signatures."
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
        "--target_signature",
        type=str,
        required=True,
        help="Target signature identifier (e.g., SBS3).",
    )
    parser.add_argument(
        "--output_plot",
        type=str,
        required=True,
        help="Path to save output heatmap plot PNG.",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=15,
        help="Number of top impacted non-target signatures to plot (default: 15).",
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="vlag",
        help="Diverging Seaborn colormap name (default: vlag).",
    )

    args = parser.parse_args()

    plot_signature_bleed_heatmap(
        file_with=args.file_with,
        file_without=args.file_without,
        target_signature=args.target_signature,
        output_plot=args.output_plot,
        top_n=args.top_n,
        cmap=args.cmap,
    )

# process PLOT_SIGNATURE_BLEED {
#     input:
#     path file_with_csv
#     path file_without_csv
#     val  target_sig
#     val  top_n_sigs

#     output:
#     path "*.png"

#     script:
#     """
#     plot_signature_bleed_heatmap.py \\
#         --file_with ${file_with_csv} \\
#         --file_without ${file_without_csv} \\
#         --target_signature ${target_sig} \\
#         --top_n ${top_n_sigs} \\
#         --output_plot "heatmap_bleed_${target_sig}.png"
#     """
# }