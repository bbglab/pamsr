import os
import pandas as pd

from .objects import MutationGroup

def load_metadata(path:str,delimiter:str):
    """
    Load and validate a metadata file.

    The function reads a metadata table from a CSV or TSV file, checks that
    the required columns are present, and verifies that there are no missing
    values in the essential fields (`sample_id` and `group`).

    Parameters
    ----------
    path : str
        Path to the metadata file to load.

    delimiter : str
        Format of the metadata file. Supported values are:
        
        - ``"csv"`` : comma-separated values (`,`)
        - ``"tsv"`` : tab-separated values (`\\t`)

    Returns
    -------
    pandas.DataFrame
        A dataframe containing the metadata information. The dataframe must
        contain at least the following columns:

        - ``sample_id`` : unique identifier for each sample.
        - ``group`` : group assignment associated with each sample.

    Raises
    ------
    ValueError
        If:
        
        - The provided delimiter type is not supported.
        - Required columns (`sample_id` or `group`) are missing.
        - Any value in the `sample_id` column is missing.
        - Any value in the `group` column is missing.

    Examples
    --------
    >>> metadata = load_metadata(
    ...     "input_data/metadata.csv",
    ...     delimiter="csv"
    ... )
    >>> metadata.head()

       sample_id group
    0   PD4199a     A
    1   PD4005a     A
    2   PD3851a     A

    Notes
    -----
    The metadata file is expected to contain one row per sample. The
    `sample_id` column is used to match samples in the mutation count matrix,
    while the `group` column is used to organize samples into independent
    analysis groups.
    """

    # Selecting the delimiter type
    delimiter_types={
        'csv':',',
        'tsv':'\t'
        }
    
    if delimiter not in delimiter_types.keys():
        raise ValueError(
            "Delimiter type is not supported"
        )
    else:
        delimiter_focus=delimiter_types[delimiter]

    metadata = pd.read_csv(path,sep=delimiter_focus)

    required_columns = [
        "sample_id",
        "group"
    ]

    missing = set(required_columns) - set(metadata.columns)

    if missing:
        raise ValueError(
            f"Metadata missing columns: {missing}"
        )


    if metadata["sample_id"].isna().any():
        raise ValueError(
            "There are missing sample_id values in metadata"
        )


    if metadata["group"].isna().any():
        raise ValueError(
            "There are missing group values in metadata"
        )


    return metadata

def load_mutation_matrix(path, delimiter):
    """
    Load and validate a mutation count matrix.

    This function reads a mutation matrix from a CSV or TSV file and performs
    basic validation checks to ensure that the matrix follows the expected
    format for SBS96 mutational signature analysis.

    The expected structure of the input file is:

    - First column: ``Mutation Types`` containing the 96 possible SBS96
      mutation channels.
    - Remaining columns: sample identifiers containing integer mutation counts.

    Parameters
    ----------
    path : str
        Path to the mutation matrix file.

    delimiter : str
        Format of the input file. Supported values are:

        - ``"csv"`` : comma-separated values (`,`)
        - ``"tsv"`` : tab-separated values (`\\t`)

    Returns
    -------
    pandas.DataFrame
        Dataframe containing the mutation count matrix. The first column
        corresponds to mutation channels and the remaining columns correspond
        to samples.

    Raises
    ------
    ValueError
        If:
        - The provided delimiter type is not supported.
        - The input file does not contain a ``Mutation Types`` column.
        - The matrix does not contain exactly 96 mutation channels.

    Notes
    -----
    The function checks that the mutation channels correspond to the SBS96
    classification scheme:

    - Six possible base substitutions:
        - C>A
        - C>G
        - C>T
        - T>A
        - T>C
        - T>G

    - Combined with the 16 possible trinucleotide contexts:

        6 substitutions x 16 contexts = 96 mutation types

    Mutation channels that do not belong to the SBS96 catalogue are reported
    but do not currently stop execution.

    Examples
    --------
    >>> mutation_matrix = load_mutation_matrix(
    ...     "input_data/in_data_a.csv",
    ...     delimiter="csv"
    ... )
    >>> mutation_matrix.head()

      Mutation Types  PD4199a  PD4005a  PD3851a
    0      A[C>A]A        58       74       31
    1      A[C>A]C        36       66       34
    2      A[C>A]G        13       12        9

    """

    delimiter_types = {
        "csv": ",",
        "tsv": "\t"
    }

    if delimiter not in delimiter_types.keys():
        raise ValueError(
            "Delimiter type is not supported. "
            "Choose between 'csv' and 'tsv'."
        )

    delimiter_focus = delimiter_types[delimiter]

    df = pd.read_csv(
        path,
        sep=delimiter_focus
    )


    if "Mutation Types" not in df.columns:
        raise ValueError(
            "Input file does not contain 'Mutation Types' column"
        )

    mutation_types = (
        'A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T',
        'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T',
        'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T',
        'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T',
        'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T',
        'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T',
        'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T',
        'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T',
        'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
        'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T',
        'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T',
        'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T',
        'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T',
        'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T',
        'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T',
        'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T',
        'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T',
        'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T',
        'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T',
        'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T',
        'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
        'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T',
        'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T',
        'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T'
    )

    if len(df.iloc[:, 0]) != 96:
        raise ValueError(
            "Input file does not contain exactly 96 mutation types"
        )

    invalid_channels = [
        val for val in df["Mutation Types"]
        if val not in mutation_types
    ]

    if invalid_channels:
        print(
            f"Found {len(invalid_channels)} channels not in the SBS96 catalogue:"
            f"\n{invalid_channels}"
        )

    return df

def load_sample_sheet(path: str, delimiter: str):
    """
    Load and validate a sample sheet file.

    The sample sheet contains the paths to the mutation matrices that will be
    analyzed and the delimiter format used by each individual mutation matrix.

    The expected structure of the sample sheet is:

    - ``input_path`` : path to a mutation matrix file.
    - ``delimiter`` : delimiter type of the corresponding mutation matrix.

    Example:

    input_path,delimiter
    input/mutational_matrices/in_data_a.csv,csv
    input/mutational_matrices/in_data_b.csv,csv
    input/mutational_matrices/in_data_c.tsv,tsv

    Parameters
    ----------
    path : str
        Path to the sample sheet file.

    delimiter : str
        Format of the sample sheet. Supported values are:

        - ``"csv"`` : comma-separated values (`,`)
        - ``"tsv"`` : tab-separated values (`\\t`)

    Returns
    -------
    pandas.DataFrame
        Dataframe containing the sample sheet information with the required
        columns:

        - ``input_path`` : path to the mutation matrix.
        - ``delimiter`` : delimiter format of the mutation matrix.

    Raises
    ------
    ValueError
        If:

        - The provided delimiter type is not supported.
        - Required columns (`input_path` or `delimiter`) are missing.

    Notes
    -----
    The delimiter specified inside the sample sheet is validated independently
    for every mutation matrix. Files containing unsupported delimiters are not
    removed from the dataframe but a warning is generated.

    This function does not load the mutation matrices themselves. The returned
    dataframe is later used by functions responsible for loading and validating
    mutation count matrices.

    Examples
    --------
    >>> sample_sheet = load_sample_sheet(
    ...     "input_data/sample_sheet.csv",
    ...     delimiter="csv"
    ... )
    >>> sample_sheet.head()

                                      input_path delimiter
    0 input/mutational_matrices/in_data_a.csv       csv
    1 input/mutational_matrices/in_data_b.csv       csv
    2 input/mutational_matrices/in_data_c.csv       csv

    """

    sample_sheet = pd.read_csv(
        path,
        sep={
            "csv": ",",
            "tsv": "\t"
        }.get(delimiter)
    )

    required_columns = [
        "input_path",
        "delimiter"
    ]

    missing = (
        set(required_columns)
        -
        set(sample_sheet.columns)
    )

    if missing:
        raise ValueError(
            f"Sample sheet missing columns: {missing}"
        )

    # Validate the delimiters specified for each mutation matrix
    valid_delimiters = {
        "csv",
        "tsv"
    }

    invalid_delimiters = sample_sheet[
        ~sample_sheet["delimiter"].isin(valid_delimiters)
    ]

    if len(invalid_delimiters) > 0:
        for _, row in invalid_delimiters.iterrows():
            print(
                "WARNING: Unsupported delimiter "
                f"{row['delimiter']} for file "
                f"{row['input_path']}"
            )

    return sample_sheet

def build_sample_registry(path):
    """
    Build a registry linking sample identifiers to their mutation matrix.

    Parameters
    ----------
    sample_sheet : pandas.DataFrame
        Dataframe containing the columns:
        - input_path
        - delimiter

    Returns
    -------
    pandas.DataFrame
        Table containing:
        - sample_id
        - matrix_path

    Notes
    -----
    Samples are extracted from the columns of every mutation matrix.
    Mutation matrices with unsupported delimiters are skipped.
    """

    valid_delimiters = {
        "csv",
        "tsv"
    }

    sample_registry = []

    for _, row in path.iterrows():
        matrix_path = row["input_path"]
        delimiter = row["delimiter"]
        if delimiter not in valid_delimiters:
            print(
                f"WARNING: Skipping {matrix_path} "
                f"because delimiter '{delimiter}' is invalid"
            )

            continue

        matrix = load_mutation_matrix(
            path=matrix_path,
            delimiter=delimiter
        )

        samples = [
            x for x in matrix.columns
            if x != "Mutation Types"
        ]

        for sample in samples:
            sample_registry.append(
                {
                    "sample_id": sample,
                    "matrix_path": matrix_path
                }
            )

    return pd.DataFrame(sample_registry)

def validate_samples(matrix, metadata):
    """
    Validate the consistency of sample identifiers between a mutation matrix
    and a metadata dataframe.

    This function checks whether all samples present in the mutation matrix
    have corresponding entries in the metadata file. Samples present in the
    metadata but absent from the mutation matrix are allowed, since metadata
    may contain information for samples that are not included in the current
    analysis.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Mutation count matrix containing sample identifiers as column names.
        The first column is expected to contain mutation channel identifiers
        (e.g., ``Mutation Types``), while the remaining columns correspond to
        samples.

    metadata : pandas.DataFrame
        Metadata dataframe containing sample information. It must contain a
        ``sample_id`` column with sample identifiers.

    Returns
    -------
    None

        The function does not return a value. Validation is performed through:

        - Raising an error when samples are present in the mutation matrix but
          absent from metadata.
        - Printing a warning when samples are present in metadata but absent
          from the mutation matrix.

    Raises
    ------
    ValueError
        If one or more samples are found in the mutation matrix that are not
        present in the metadata dataframe.

    Notes
    -----
    The validation rules are asymmetric:

    - Mutation matrix → metadata:
        All samples must have metadata information because downstream analyses
        require group assignments and sample annotations.

    - Metadata → mutation matrix:
        Missing samples are allowed because metadata files may contain samples
        that are not included in the current dataset.

    Examples
    --------
    >>> validate_samples(
    ...     matrix,
    ...     metadata
    ... )

    If the mutation matrix contains:

    ``PD4199a, PD4005a, PD9999a``

    and metadata contains:

    ``PD4199a, PD4005a``

    the function raises:

    ValueError:
        Samples in input but not metadata: {'PD9999a'}

    If metadata contains additional samples:

    ``PD4199a, PD4005a, PD8888a``

    and the matrix contains:

    ``PD4199a, PD4005a``

    a warning is printed, but execution continues.

    """

    input_samples = set(matrix.columns[1:])
    metadata_samples = set(metadata["sample_id"])

    # Samples in mutation matrix but not metadata
    missing_metadata = input_samples - metadata_samples

    if missing_metadata:
        raise ValueError(
            f"Samples in input but not metadata: {missing_metadata}"
        )

    # Samples in metadata but not input are allowed
    missing_input = metadata_samples - input_samples

    if missing_input:
        print(
            f"Warning: Samples in metadata but not input: {missing_input}"
        )

def validate_duplicate_samples(sample_registry):
    """
    Detect duplicated samples appearing in multiple mutation matrices.

    Parameters
    ----------
    sample_registry : pandas.DataFrame
        Registry containing:
        
        - ``sample_id`` : sample identifier
        - ``matrix_path`` : path to the mutation matrix containing the sample

    Returns
    -------
    pandas.DataFrame
        Registry after removing duplicated samples.

    Notes
    -----
    If a sample appears in multiple mutation matrices, all occurrences are
    removed from the registry to avoid ambiguous sample assignment.

    The function reports:
    
    - duplicated sample identifier
    - number of mutation matrices where it appears
    - paths of all mutation matrices containing that sample

    Examples
    --------
    If:

    sample_id | matrix_path
    ----------|------------
    PD4199a   | in_data_a.csv
    PD4199a   | in_data_b.csv
    PD4199a   | in_data_c.csv

    the warning will indicate that PD4199a appears in 3 datasets.
    """

    duplicated_samples = (
        sample_registry
        .groupby("sample_id")
        .filter(lambda x: len(x) > 1)
    )

    if not duplicated_samples.empty:

        print(
            "\nWARNING: duplicated samples found:"
        )

        duplicated_groups = (
            duplicated_samples
            .groupby("sample_id")
        )

        for sample_id, locations in duplicated_groups:

            datasets = locations["matrix_path"].tolist()

            print(
                f"\nSample: {sample_id}"
            )

            print(
                f"Found in {len(datasets)} mutation matrices:"
            )

            for dataset in datasets:
                print(
                    f"  - {dataset}"
                )

        # Remove duplicated samples completely
        sample_registry = (
            sample_registry
            .drop_duplicates(
                subset="sample_id",
                keep=False
            )
        )

    return sample_registry

def validate_registry_against_metadata(sample_registry, metadata):
    """
    Validate the consistency between a sample registry and the metadata file.

    This function compares the sample identifiers present in the sample
    registry against those listed in the metadata file.

    Two independent validations are performed:

    1. Samples present in one or more mutation matrices but absent from the
       metadata file. These samples are reported together with the mutation
       matrix in which they were found.

    2. Samples present in the metadata file but absent from all mutation
       matrices. These samples are reported as warnings since metadata may
       intentionally contain samples that are not included in the current
       analysis.

    Parameters
    ----------
    sample_registry : pandas.DataFrame
        Dataframe linking every sample to the mutation matrix where it was
        found. The dataframe must contain at least the following columns:

        - ``sample_id`` : unique sample identifier.
        - ``matrix_path`` : path to the mutation matrix containing the sample.

    metadata : pandas.DataFrame
        Metadata dataframe containing information for every sample. The
        dataframe must contain at least the column:

        - ``sample_id`` : unique sample identifier.

    Returns
    -------
    None
        The function does not return any value.

        Validation results are communicated through printed warnings.

    Notes
    -----
    The validation is asymmetric.

    Samples present in mutation matrices but absent from the metadata file
    usually indicate an error in the input files and should be corrected before
    downstream analyses.

    Samples present in the metadata but absent from mutation matrices are
    allowed. This situation commonly occurs when the metadata describes more
    samples than those included in the current analysis.

    The function assumes that duplicated samples have already been removed by
    ``validate_duplicate_samples()``.

    Examples
    --------
    >>> validate_registry_against_metadata(
    ...     sample_registry,
    ...     metadata
    ... )

    Example output:

    WARNING: Samples present in mutation matrices but missing from metadata:

        PD4199a -> input/mutational_matrices/in_data_a.csv
        PD4005a -> input/mutational_matrices/in_data_b.csv

    WARNING: Samples present in metadata but absent from mutation matrices:

        PD5001a
        PD5002a
    """

    registry_samples = set(
        sample_registry["sample_id"]
    )

    metadata_samples = set(
        metadata["sample_id"]
    )

    # Samples present in mutation matrices but absent from metadata

    missing_metadata = (
        registry_samples -
        metadata_samples
    )

    if missing_metadata:

        print(
            "\nWARNING: Samples present in mutation "
            "matrices but missing from metadata:"
        )

        for sample in missing_metadata:

            location = sample_registry.loc[
                sample_registry["sample_id"] == sample,
                "matrix_path"
            ].values[0]

            print(
                f"  {sample} -> {location}"
            )

    # Samples present in metadata but absent from mutation matrices
    missing_matrix = (
        metadata_samples -
        registry_samples
    )

    if missing_matrix:
        print(
            "\nWARNING: Samples present in metadata "
            "but absent from mutation matrices:"
        )
        for sample in sorted(missing_matrix):

            print(
                f"  {sample}"
            )

def create_groups(sample_sheet, metadata):
    """
    Create one MutationGroup object per metadata group.

    Returns
    -------
    dict
        Dictionary indexed by group identifier.
    """

    groups = {}

    for _, row in sample_sheet.iterrows():

        matrix = load_mutation_matrix(
            row["input_path"],
            row["delimiter"]
        )

        for group_id, group_metadata in metadata.groupby("group"):

            sample_ids = [
                s
                for s in group_metadata["sample_id"]
                if s in matrix.columns
            ]

            if len(sample_ids) == 0:
                continue

            group_matrix = matrix[
                sample_ids
            ]

            if group_id not in groups:

                groups[group_id] = MutationGroup(
                    group_id,
                    group_matrix,
                    group_metadata,
                    matrix["Mutation Types"]
                )

            else:

                groups[group_id].samples = (
                    groups[group_id]
                    .samples
                    .join(group_matrix)
                )

    return groups

def create_groups_from_sample_sheet(sample_sheet):
    """
    Create one MutationGroup object per group defined directly in the sample sheet.

    This function bypasses the need for an external metadata file by grouping
    samples based on a 'group' or 'model' column provided directly in the 
    sample sheet.

    Parameters
    ----------
    sample_sheet : pandas.DataFrame
        Dataframe containing the sample sheet information. It must contain:
        - `input_path`
        - `delimiter`
        - `group` (or `model`)

    Returns
    -------
    dict
        Dictionary of MutationGroup objects indexed by group identifier.
    """
    groups = {}

    # Identify the column used for grouping in the sample sheet
    if "group" in sample_sheet.columns:
        group_col = "group"
    elif "model" in sample_sheet.columns:
        group_col = "model"
    else:
        raise ValueError(
            "Sample sheet must contain a 'group' or 'model' column to define groups."
        )

    for _, row in sample_sheet.iterrows():
        # Load the mutation matrix for the current sample sheet row
        matrix = load_mutation_matrix(
            row["input_path"],
            row["delimiter"]
        )

        group_id = row[group_col]

        # Extract sample columns (ignoring the 'Mutation Types' column)
        sample_ids = [
            s for s in matrix.columns 
            if s != "Mutation Types"
        ]

        if len(sample_ids) == 0:
            continue

        group_matrix = matrix[sample_ids]

        # Create a mock metadata dataframe to satisfy MutationGroup requirements
        mock_metadata = pd.DataFrame({
            "sample_id": sample_ids,
            "group": group_id
        })

        if group_id not in groups:
            groups[group_id] = MutationGroup(
                group_id,
                group_matrix,
                mock_metadata,
                matrix["Mutation Types"]
            )
        else:
            # Append new samples to the existing group matrix
            groups[group_id].samples = (
                groups[group_id]
                .samples
                .join(group_matrix)
            )

            # Append to mock metadata
            groups[group_id].metadata = pd.concat(
                [groups[group_id].metadata, mock_metadata], 
                ignore_index=True
            )

    return groups