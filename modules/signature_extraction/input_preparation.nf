process INPUT_PREPARATION {
    tag "$metadata.name"
    // conda = "/home/alberto/miniforge3/envs/gsd"

    // Mount local bin directory into the container's PYTHONPATH
    // env PYTHONPATH: "${projectDir}/bin"

    publishDir "${params.outdir}/signature_extraction/input_prep", mode: 'copy'

    input:
    path metadata
    path samplesheet
    val  metadata_delim
    val  samplesheet_delim

    output:
    path "validated_registry.csv", emit: registry

    script:
    """
    #!/usr/bin/env python3
    import os

    # Link project input directory into the process workspace for relative path resolution
    if not os.path.exists("input"):
        os.symlink("${projectDir}/input", "input")

    from help_utils.data_loader_validator import (
        load_metadata,
        load_sample_sheet,
        build_sample_registry,
        validate_duplicate_samples,
        validate_registry_against_metadata,
        create_groups
    )

    metadata_df = load_metadata("${metadata}", "${metadata_delim}")
    sample_sheet_df = load_sample_sheet("${samplesheet}", "${samplesheet_delim}")

    registry = build_sample_registry(sample_sheet_df)
    registry = validate_duplicate_samples(registry)

    validate_registry_against_metadata(registry, metadata_df)
    groups = create_groups(sample_sheet_df, metadata_df)

    registry.to_csv("validated_registry.csv", index=False)
    """
}