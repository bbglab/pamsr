#!/usr/bin/env python3
"""
nnls_signature_check.py

Check how well a set of "background" mutational signatures can explain a
"target" signature, using non-negative least squares (NNLS).

This is useful *before* running a costly signature-detection power analysis:
if a target signature can already be closely approximated by a non-negative
mixture of the background signatures you plan to use in a presence test
(e.g. mSigAct's SignaturePresenceTest), you should expect reduced /
non-monotonic detection power at low-to-moderate injected exposure levels,
because the "without target" null model can partially reallocate exposure
among the background signatures to mimic the target.

Input
-----
A signature catalog file (.tsv or .csv), with:
  - first column: mutation channel / type labels (e.g. 96 SBS trinucleotide
    contexts), used as the row index
  - remaining columns: one column per signature, values proportional to
    probability (need not already sum to 1 -- they are renormalized)

Usage
-----
    python3 nnls_signature_check.py \\
        --catalog COSMIC_v3_6_SBS_GRCh38.txt \\
        --target SBS31 \\
        --background SBS1,SBS5,SBS40a

    # Explicit delimiter override (otherwise inferred from file extension,
    # falling back to auto-sniffing) and saving the residual profile:
    python3 nnls_signature_check.py \\
        --catalog my_catalog.csv \\
        --delimiter , \\
        --target SBS31 \\
        --background SBS1,SBS5,SBS40a \\
        --output-residual residual_SBS31.csv

Output
------
Printed to stdout:
  - NNLS mixture coefficients (raw and normalized to sum to 1)
  - Fraction of the target signature's probability mass explainable by the
    non-negative background mixture (1 - total variation distance)
  - L1 / total-variation distance and cosine similarity between the target
    and the best-fit background mixture
  - Pairwise cosine similarity between the target and each individual
    background signature (helps identify which one(s) drive collinearity)

Optionally, --output-residual writes a CSV with per-channel target,
fitted, and residual values.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import nnls


def infer_delimiter(path: str, delimiter_arg: str | None) -> str:
    """Decide which delimiter to use for reading the catalog file."""
    if delimiter_arg:
        # Allow convenient aliases
        aliases = {"tab": "\t", "\\t": "\t", "comma": ",", ",": ",", "tsv": "\t", "csv": ","}
        return aliases.get(delimiter_arg.lower(), delimiter_arg)

    ext = os.path.splitext(path)[1].lower()
    if ext == ".tsv":
        return "\t"
    if ext == ".csv":
        return ","
    # Fall back: sniff the first line
    with open(path, "r") as fh:
        first_line = fh.readline()
    if "\t" in first_line:
        return "\t"
    return ","


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return float("nan")
    return float(np.dot(a, b) / denom)


def load_catalog(path: str, delimiter: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=delimiter, index_col=0)
    # Be forgiving about accidental whitespace in column names
    df.columns = [str(c).strip() for c in df.columns]
    return df


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    sums = df.sum(axis=0)
    zero_sum_cols = sums[sums <= 0].index.tolist()
    if zero_sum_cols:
        raise ValueError(
            f"The following columns sum to zero and cannot be normalized: {zero_sum_cols}"
        )
    return df.divide(sums, axis=1)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run NNLS of a target signature against a set of background "
            "signatures extracted from columns of a catalog file (tsv/csv)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--catalog", required=True,
        help="Path to the signature catalog file (.tsv or .csv). "
             "First column = mutation channel labels, other columns = signatures.",
    )
    parser.add_argument(
        "--target", required=True,
        help="Name of the column to use as the target signature (e.g. SBS31).",
    )
    parser.add_argument(
        "--background", required=True,
        help="Comma-separated list of column names to use as background "
             "signatures (e.g. SBS1,SBS5,SBS40a). The target signature does "
             "NOT need to be included here.",
    )
    parser.add_argument(
        "--delimiter", default=None,
        help="Field delimiter to use when reading the catalog. If omitted, "
             "it is inferred from the file extension (.tsv -> tab, .csv -> "
             "comma), falling back to sniffing the first line. Accepts "
             "'tab'/'comma' as convenient aliases.",
    )
    parser.add_argument(
        "--output-residual", default=None,
        help="Optional path to write a CSV with per-channel target, fitted, "
             "and residual values.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.catalog):
        sys.exit(f"ERROR: catalog file not found: {args.catalog}")

    delimiter = infer_delimiter(args.catalog, args.delimiter)
    df = load_catalog(args.catalog, delimiter)

    background_names = [s.strip() for s in args.background.split(",") if s.strip()]
    target_name = args.target.strip()

    missing = [s for s in background_names + [target_name] if s not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: the following signature(s) were not found as columns in "
            f"'{args.catalog}': {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    if target_name in background_names:
        print(
            f"NOTE: target signature '{target_name}' is also listed in "
            f"--background; it will be excluded from the background matrix "
            f"used to explain itself.",
            file=sys.stderr,
        )
        background_names = [s for s in background_names if s != target_name]
        if not background_names:
            sys.exit("ERROR: no background signatures left after excluding the target.")

    df_norm = normalize_columns(df[[target_name] + background_names].astype(float))

    t = df_norm[target_name].to_numpy()
    B = df_norm[background_names].to_numpy()

    coef, _ = nnls(B, t)
    fitted = B @ coef
    residual = t - fitted

    l1_residual = float(np.abs(residual).sum())
    total_variation = l1_residual / 2.0
    explainable_fraction = 1.0 - total_variation
    cos_target_fitted = cosine_similarity(t, fitted)

    coef_sum = coef.sum()
    coef_normalized = coef / coef_sum if coef_sum > 0 else coef

    print("=" * 70)
    print(f"Target signature:      {target_name}")
    print(f"Background signatures: {background_names}")
    print(f"Catalog:               {args.catalog}  (delimiter={delimiter!r})")
    print(f"Number of channels:    {len(t)}")
    print("=" * 70)

    print("\nNNLS mixture coefficients (raw / normalized to sum to 1):")
    for name, c_raw, c_norm in zip(background_names, coef, coef_normalized):
        print(f"  {name:<12s} raw = {c_raw:8.4f}   normalized = {c_norm:8.4f}")

    print(f"\nTotal variation distance (target vs. best-fit background mixture): {total_variation:.4f}")
    print(f"Fraction of target explainable by background mixture:             {explainable_fraction:.4f}")
    print(f"Cosine similarity (target vs. best-fit background mixture):       {cos_target_fitted:.4f}")

    print("\nPairwise cosine similarity (target vs. each individual background signature):")
    for name in background_names:
        b = df_norm[name].to_numpy()
        print(f"  {target_name} vs {name:<12s} = {cosine_similarity(t, b):.4f}")

    if explainable_fraction > 0.85:
        interpretation = (
            "HIGH collinearity: the background set can nearly reproduce the "
            "target on its own. Expect weak / unstable / non-monotonic "
            "detection power at low-to-moderate injected exposure in a "
            "presence-test power analysis."
        )
    elif explainable_fraction > 0.5:
        interpretation = (
            "MODERATE collinearity: expect a real transition / dip zone in "
            "the detection curve at low-to-moderate injected exposure, "
            "with reliable detection only once injected exposure produces "
            "enough of the non-explainable ('distinctive') signal."
        )
    else:
        interpretation = (
            "LOW collinearity: the target signature is largely distinguishable "
            "from the background set; detection power should increase fairly "
            "smoothly with injected exposure."
        )
    print(f"\nInterpretation: {interpretation}")

    if args.output_residual:
        out_df = pd.DataFrame(
            {
                "channel": df_norm.index,
                "target": t,
                "fitted_from_background": fitted,
                "residual": residual,
            }
        )
        out_df.to_csv(args.output_residual, index=False)
        print(f"\nPer-channel residual profile written to: {args.output_residual}")


if __name__ == "__main__":
    main()