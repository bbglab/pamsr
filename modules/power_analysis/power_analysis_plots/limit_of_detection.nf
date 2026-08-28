process PLOT_LIMIT_OF_DETECTION {
    input:
    path stats_csv
    val  alpha_threshold

    output:
    path "*.png"

    script:
    """
    plot_limit_of_detection.py \\
        --csv_file ${stats_csv} \\
        --alpha ${alpha_threshold} \\
        --output_plot "limit_of_detection_p${alpha_threshold}.png"
    """
}