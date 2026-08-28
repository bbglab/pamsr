process COMPARE_EXPOSURES {
    input:
    path file_with_csv
    path file_without_csv
    val  target_signature

    output:
    path "*.png"

    script:
    """
    compare_exposures.py \\
        --file_with ${file_with_csv} \\
        --file_without ${file_without_csv} \\
        --target_sig ${target_signature} \\
        --output_plot "comparison_${target_signature}.png"
    """
}