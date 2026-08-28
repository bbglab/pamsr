#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression


def plot_limit_of_detection(
    csv_file: str,
    output_plot: str,
    alpha: float = 0.05,
    color: str = "crimson",
) -> None:
    """
    Fits a Logistic Regression model to derive Limit of Detection (LoD50 and LoD95) 
    curves for signature presence based on chi-square p-values.

    Parameters
    ----------
    csv_file : str
        Path to input CSV containing `n_injected` and `chisq_p`.
    output_plot : str
        File path where the output PNG plot will be saved.
    alpha : float, default=0.05
        Significance threshold for binary detection classification (`chisq_p < alpha`).
    color : str, default="crimson"
        Color for fitted logistic sigmoid curve.

    Returns
    -------
    None
        Saves output figure directly to `output_plot`.

    Raises
    ------
    FileNotFoundError
        If `csv_file` does not exist at the specified path.
    KeyError
        If mandatory columns `n_injected` or `chisq_p` are missing.

    Notes
    -----
    - Uses logit transformation $X = \\frac{\\text{logit}(p) - b_0}{b_1}$ to solve for $LoD_{50}$ ($p=0.50$) and $LoD_{95}$ ($p=0.95$).
    - Overlays empirical detection proportions per injection level alongside fitted curve.

    Examples
    --------
    >>> plot_limit_of_detection(
    ...     csv_file="output/msigact/summary_stats.csv",
    ...     output_plot="output/msigact/limit_of_detection_5000.png",
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
    df = df.dropna(subset=["n_injected", "chisq_p"])

    # Binary outcome: 1 if significant target detection, 0 otherwise
    X = df[["n_injected"]].values
    y = (df["chisq_p"] < alpha).astype(int).values

    # Fit Logistic Regression model
    model = LogisticRegression()
    model.fit(X, y)

    b0 = model.intercept_[0]
    b1 = model.coef_[0][0]

    # Calculate LoD50 and LoD95 using logit formula: X = (logit(p) - b0) / b1
    lod_50 = (np.log(0.50 / (1 - 0.50)) - b0) / b1 if b1 > 0 else np.nan
    lod_95 = (np.log(0.95 / (1 - 0.95)) - b0) / b1 if b1 > 0 else np.nan

    # Generate fine grid for smooth sigmoid curve
    x_range = np.linspace(X.min(), X.max(), 500).reshape(-1, 1)
    y_prob = model.predict_proba(x_range)[:, 1]

    # Compute empirical detection rates for scatter overlay
    empirical = (
        df.groupby("n_injected")["chisq_p"]
        .apply(lambda p: (p < alpha).mean())
        .reset_index()
    )

    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid")

    # Empirical data points
    plt.scatter(
        empirical["n_injected"],
        empirical["chisq_p"],
        color="black",
        s=40,
        zorder=3,
        label="Empirical Data",
    )

    # Fitted Logistic Sigmoid Curve
    plt.plot(
        x_range,
        y_prob,
        color=color,
        linewidth=2.5,
        label="Logistic Fit (Sigmoid)",
    )

    # Threshold markers
    if not np.isnan(lod_50) and lod_50 >= X.min():
        plt.axvline(
            lod_50,
            color="goldenrod",
            linestyle="--",
            linewidth=1.5,
            label=f"$LoD_{{50}}$: {lod_50:.1f} muts",
        )
        plt.axhline(0.50, color="goldenrod", linestyle=":", alpha=0.6)

    if not np.isnan(lod_95) and lod_95 <= X.max() * 1.2:
        plt.axvline(
            lod_95,
            color="forestgreen",
            linestyle="--",
            linewidth=1.5,
            label=f"$LoD_{{95}}$: {lod_95:.1f} muts",
        )
        plt.axhline(0.95, color="forestgreen", linestyle=":", alpha=0.6)

    plt.title(
        f"Limit of Detection ($LoD$) Curve (Significance Threshold $\\alpha = {alpha}$)",
        fontsize=13,
    )
    plt.xlabel(r"Number of Injected Mutations ($n\_injected$)", fontsize=11)
    plt.ylabel(r"Detection Probability ($P(p < \alpha)$)", fontsize=11)
    plt.ylim(-0.05, 1.05)
    plt.legend(loc="lower right", fontsize=10)
    plt.tight_layout()

    # Save output plot
    out_dir = os.path.dirname(output_plot)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit logistic regression to estimate LoD50 and LoD95 signature detection limits."
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
        help="Significance threshold alpha (default: 0.05).",
    )
    parser.add_argument(
        "--color",
        type=str,
        default="crimson",
        help="Colormap curve color (default: crimson).",
    )

    args = parser.parse_args()

    plot_limit_of_detection(
        csv_file=args.csv_file,
        output_plot=args.output_plot,
        alpha=args.alpha,
        color=args.color,
    )

    

# process PLOT_LIMIT_OF_DETECTION {
#     input:
#     path stats_csv
#     val  alpha_threshold

#     output:
#     path "*.png"

#     script:
#     """
#     plot_limit_of_detection.py \\
#         --csv_file ${stats_csv} \\
#         --alpha ${alpha_threshold} \\
#         --output_plot "limit_of_detection_p${alpha_threshold}.png"
#     """
# }