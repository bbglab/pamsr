process PLOT_RECONSTRUCTION_RESIDUALS {
    input:
    path stats_csv

    output:
    path "*.png"

    script:
    """
    plot_reconstruction_residuals.py \\
        --csv_file ${stats_csv} \\
        --output_plot "reconstruction_residuals.png"
    """
}