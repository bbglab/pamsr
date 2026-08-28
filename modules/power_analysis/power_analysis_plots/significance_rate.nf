process SIGNIFICANCE_RATE {
    input:
    path stats_csv
    val  alpha_val

    output:
    path "*.png"

    script:
    """
    plot_significance_rate.py \\
        --csv_file ${stats_csv} \\
        --alpha ${alpha_val} \\
        --output_plot "significance_rate_p${alpha_val}.png"
    """
}
