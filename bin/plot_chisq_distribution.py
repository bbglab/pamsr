#!/usr/bin/env python3
"""
plot_chisq_distribution.py

Diagnose calibration of a SignaturePresenceTest-style likelihood-ratio test
across an injected-mutation grid, using `merged_summary_statistics.csv` from
the PAMSR power-analysis output (columns expected: n_injected, statistic,
chisq_p -- extra columns are ignored).

Produces three figures:

1. qqplot_null_pvalues.png
   QQ-plot of -log10(p) at the "null" injection level (default: n_injected=0)
   against the -log10(p) expected under Uniform(0,1). If the test is well
   calibrated, points should fall on the y = x diagonal. Points bowing far
   above the diagonal indicate the null is anti-conservative (e.g. because
   the background noise model is too tight, or the "without target" model
   can't fit the true baseline) -- exactly what produces "100% significant
   at every injection level" symptoms.

2. statistic_histograms.png
   Small-multiple histograms of the LRT `statistic`, one panel per selected
   injection level, overlaid with:
     - the theoretical chi-square(df=1) density (textbook LRT reference)
     - the 50:50 mixture of chi-square(df=0) [point mass at 0] and
       chi-square(df=1) (the "chi-bar-square" reference distribution that
       applies when the null parameter sits on a boundary, e.g. a
       non-negativity constraint on signature exposure -- this is usually
       the more appropriate reference for signature attribution LRTs).

3. statistic_boxplot_by_injection.png
   Boxplot/violin of the `statistic` distribution across the full injected
   mutation grid, to see the whole spread evolve (not just the % passing
   a threshold).

Usage
-----
    python3 plot_chisq_distribution.py \\
        --summary_csv merged_summary_statistics.csv \\
        --output_dir plots/ \\
        --null_level 0 \\
        --levels 0,100,300,600,1000

If --levels is omitted, a spread of up to 6 levels across the grid
(including the null level) is chosen automatically.
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def load_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"n_injected", "statistic", "chisq_p"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required column(s) in {path}: {missing}")
    df["n_injected"] = pd.to_numeric(df["n_injected"], errors="coerce")
    df["statistic"] = pd.to_numeric(df["statistic"], errors="coerce")
    df["chisq_p"] = pd.to_numeric(df["chisq_p"], errors="coerce")
    df = df.dropna(subset=["n_injected", "statistic", "chisq_p"])
    return df


def pick_levels(all_levels, null_level, n_panels=6):
    all_levels = sorted(all_levels)
    if null_level not in all_levels:
        # still show it if it's close, otherwise just proceed without it
        pass
    if len(all_levels) <= n_panels:
        return all_levels
    # evenly spaced indices, always keep the null level and the max level
    idx = np.linspace(0, len(all_levels) - 1, n_panels).round().astype(int)
    chosen = sorted(set(all_levels[i] for i in idx))
    if null_level in all_levels and null_level not in chosen:
        chosen = sorted(chosen + [null_level])
    return chosen


def qqplot_null_pvalues(df: pd.DataFrame, null_level: float, output_path: str):
    sub = df.loc[df["n_injected"] == null_level, "chisq_p"].to_numpy()
    if len(sub) == 0:
        print(f"WARNING: no rows found at n_injected == {null_level}; skipping QQ-plot.")
        return

    n = len(sub)
    observed = np.sort(sub)
    # avoid log(0)
    eps = 1e-300
    observed = np.clip(observed, eps, 1.0)
    expected = (np.arange(1, n + 1) - 0.5) / n

    obs_log = -np.log10(observed)
    exp_log = -np.log10(expected)

    # Kolmogorov-Smirnov test against Uniform(0,1) as a numeric calibration check
    ks_stat, ks_p = stats.kstest(sub, "uniform")
    frac_below_05 = float((sub < 0.05).mean())

    plt.figure(figsize=(6, 6))
    max_val = max(obs_log.max(), exp_log.max()) * 1.05
    plt.plot([0, max_val], [0, max_val], color="grey", linestyle="--", linewidth=1, label="expected (y = x)")
    plt.scatter(exp_log, obs_log, s=14, alpha=0.6, color="crimson")
    plt.xlabel("Expected -log10(p) under Uniform(0,1)")
    plt.ylabel("Observed -log10(p)")
    plt.title(
        f"QQ-plot of null p-values (n_injected = {null_level}, n = {n})\n"
        f"KS test vs Uniform(0,1): stat = {ks_stat:.3f}, p = {ks_p:.2e}\n"
        f"Fraction with p < 0.05: {frac_below_05:.3f} (expect ~0.05 if calibrated)"
    )
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")
    print(
        f"  Null calibration check (n_injected={null_level}): "
        f"KS p={ks_p:.2e}, fraction p<0.05 = {frac_below_05:.3f} "
        f"(should be close to the alpha you use, ~0.05, if calibrated)"
    )


def statistic_histograms(df: pd.DataFrame, levels, output_path: str):
    n_panels = len(levels)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.8 * nrows), squeeze=False)

    x = np.linspace(0.001, 30, 500)
    chi2_pdf = stats.chi2.pdf(x, df=1)

    for i, level in enumerate(levels):
        ax = axes[i // ncols][i % ncols]
        sub = df.loc[df["n_injected"] == level, "statistic"].to_numpy()
        sub = sub[np.isfinite(sub)]
        if len(sub) == 0:
            ax.set_title(f"n_injected = {level} (no data)")
            continue

        x_max = max(10, np.percentile(sub, 99) * 1.2)
        x_local = np.linspace(0.001, x_max, 500)
        chi2_local = stats.chi2.pdf(x_local, df=1)
        # 50:50 mixture of chi2(df=0) [point mass at 0] and chi2(df=1):
        # represented here as 0.5 * chi2(df=1) density for the continuous part,
        # with the remaining 0.5 mass sitting exactly at statistic = 0.
        mixture_continuous = 0.5 * chi2_local

        ax.hist(sub, bins=40, density=True, alpha=0.55, color="steelblue", label="observed")
        ax.plot(x_local, chi2_local, color="crimson", linewidth=1.5, label="chi-sq(df=1)")
        ax.plot(x_local, mixture_continuous, color="darkorange", linewidth=1.5,
                linestyle="--", label="0.5*chi-sq(df=1) (+0.5 mass at 0)")
        ax.set_title(f"n_injected = {level}  (n = {len(sub)})")
        ax.set_xlabel("LRT statistic")
        ax.set_ylabel("density")
        ax.legend(fontsize=7)

    # hide unused panels
    for j in range(n_panels, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")


def statistic_boxplot(df: pd.DataFrame, output_path: str):
    levels = sorted(df["n_injected"].unique())
    data = [df.loc[df["n_injected"] == lv, "statistic"].to_numpy() for lv in levels]

    plt.figure(figsize=(max(8, 0.5 * len(levels)), 5))
    plt.boxplot(data, positions=range(len(levels)), showfliers=False, widths=0.6)
    plt.xticks(range(len(levels)), [str(int(lv)) if float(lv).is_integer() else str(lv) for lv in levels],
               rotation=45, ha="right")
    plt.axhline(stats.chi2.ppf(0.95, df=1), color="grey", linestyle="--", linewidth=1,
                label="chi-sq(1) 0.95 quantile (~3.84)")
    plt.xlabel("n_injected")
    plt.ylabel("LRT statistic")
    plt.title("Distribution of the LRT statistic across the injection grid")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary_csv", required=True,
                         help="Path to merged_summary_statistics.csv (needs n_injected, statistic, chisq_p).")
    parser.add_argument("--output_dir", default="chisq_diagnostics",
                         help="Directory to write the output figures (created if missing).")
    parser.add_argument("--null_level", type=float, default=0.0,
                         help="Value of n_injected to treat as the signal-free null (default: 0).")
    parser.add_argument("--levels", default=None,
                         help="Comma-separated list of n_injected values to show in the histogram panel. "
                              "If omitted, up to 6 levels are chosen automatically across the grid.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = load_summary(args.summary_csv)

    all_levels = sorted(df["n_injected"].unique())
    if args.levels:
        levels = [float(x) for x in args.levels.split(",")]
    else:
        levels = pick_levels(all_levels, args.null_level)

    print(f"Loaded {len(df)} rows across {len(all_levels)} injection levels: {all_levels}")
    print(f"Using levels for histogram panel: {levels}\n")

    qqplot_null_pvalues(df, args.null_level, os.path.join(args.output_dir, "qqplot_null_pvalues.png"))
    statistic_histograms(df, levels, os.path.join(args.output_dir, "statistic_histograms.png"))
    statistic_boxplot(df, os.path.join(args.output_dir, "statistic_boxplot_by_injection.png"))


if __name__ == "__main__":
    main()